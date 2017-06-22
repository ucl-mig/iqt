function compute_dti_respairs(input_dir, output_dir, data_folders, sub_path, ...
                     ds_rate, dw_file, bvals_file, bvecs_file, ...
                     mask_file, grad_file, dt_pref)
% COMPUTE_DTI_RESPAIRS
%   Computes DTIs on the original DWIs and its downsampled version.
%   As a result, we obtain high-res and low-res DTIs.
%
%   Args: 
%       INPUT_DIR: Input root folder (eg HCP root)
%       OUTPUT_DIR: Output root folder
%       DATA_FOLDERS: Subjects to process (eg HCP subjects)
%       SUB_PATH: Input/Output data access is done as:
%                  input_dir/data_folder/sub_path/file_name
%                  output_dir/data_folder/sub_path/file_name
%       DS_RATE: Super-resolution factor
%       DW_FILE: DWI data file
%       BVALS_FILE: b-values file
%       BVECS_FILE: b-vectors file
%       MASK_FILE: mask file
%       GRAD_FILE: grad file for HCP data (otherwise = '')
%       DT_PREF: DTI output prefix
%
%   (always end directory paths with a forward/back slash)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

check_path(input_dir);
check_path(output_dir);
check_path(sub_path);

parfor fi = 1:length(data_folders)
    if(~exist([output_dir data_folders{fi}  '/' sub_path], 'dir'))
        mkdir(output_dir, [data_folders{fi}  '/' sub_path]);
    end
    
    % -------- Step 1: Compute DTI for the original high-res DWI ---------
    fprintf('DTI computation subject: %s\n', data_folders{fi});
    input_folder = [input_dir data_folders{fi} '/' sub_path];
    output_folder = [output_dir data_folders{fi} '/' sub_path];
    
    % Read in the DWI data
    fprintf('Loading DWI: %s (%s)\n', dw_file, data_folders{fi});
    tmp = load_nii( [input_folder dw_file] );
    dw = flipdim(tmp.img,1);
    [XSIZE,YSIZE,ZSIZE,~] = size(dw);
    hdr = tmp.hdr;
    
    % Read bvals and bvecs text files
    fprintf('Loading bvals/bvecs: %s, %s (%s)\n', bvals_file, bvecs_file, data_folders{fi});
    bvecs = load([input_folder bvecs_file]); % should be 3xN
    bvals = load([input_folder bvals_file]); % should be 1xN

    % Read the brain mask
    fprintf('Loading mask: %s (%s)\n', mask_file, data_folders{fi});
    tmp = load_nii( [input_folder '/' mask_file] );
    mask = flipdim(tmp.img,1);

    % Read gradient nonlinearity file and qform (HCP only).
    g = [];
    if ~strcmp(grad_file, '')
        fprintf('Loading grad-nonlinearity: %s (%s)\n', grad_file, data_folders{fi});
        tmp = load_nii( [input_folder '/' grad_file] );
        g = flipdim(tmp.img,1);
    end
    tmp = []; % clear to run within parfor.

    % Here we want to fit only to the inner, b=1000, shell as well as b=0
    % images.
    dw_inds = find(bvals<1500);
    bvals = bvals(dw_inds);
    bvecs = bvecs(:,dw_inds);
    dw = dw(:,:,:,dw_inds);

    % Identify b=0 indices
    b0_thresh = 100;
    b0_inds = find(bvals<b0_thresh);

    % Compute DTIs.
    % Loop over and reconstruct.
    dt = zeros(XSIZE,YSIZE,ZSIZE,8);
    dt(:,:,:,1) = -1;
    dt_mask_invalid_tensors = zeros(XSIZE,YSIZE,ZSIZE);

    % Note that each voxel stores eight numbers - two more than the six
    % diffusion tensor elements.  The first entry is a code indicating
    % background (-1) or foreground (0).  The second is the log of the b=0
    % signal.  The rest are the elements of the diffusion tensor ordered as
    % (Dxx Dxy Dxz Dyy Dyz Dzz) where the actual tensor is
    % Dxx Dxy Dxz
    % Dxy Dyy Dyz
    % Dxz Dyz Dzz

    for k=1:ZSIZE
        display(sprintf('Slice %i of %i (%s).', k, ZSIZE, data_folders{fi}));
        for j=1:YSIZE
            for i=1:XSIZE
                if(mask(i,j,k))    
                    new_bvecs = bvecs;
                    new_bvals = bvals;
% %                     % Average the gradient correction over the block
% %                     L = reshape(squeeze(mean(mean(mean(g(i:(i+ds_rate-1),...
% %                                                          j:(j+ds_rate-1),...
% %                                                          k:(k+ds_rate-1),:)...
% %                                                          )))),3,3);
% %                     I = eye(3);
% % 
% %                     % correct bvecs and calculate their norm
% %                     v = (I+L)*bvecs;
% %                     n = sqrt(sum(v.^2));
% % 
% %                     % rotate bvecs to voxel coord system 
% %                     % and correct both bvecs and bvals
% %                     new_bvecs = v./repmat(n,3,1);
% %                     new_bvals = n.^2.*bvals;

                    % Construct design matrix for linear estimation of DT.
                    X = [ones(1, length(new_bvals));...
                        -new_bvals.*new_bvecs(1,:).*new_bvecs(1,:);...
                        -2*new_bvals.*new_bvecs(1,:).*new_bvecs(2,:);...
                        -2*new_bvals.*new_bvecs(1,:).*new_bvecs(3,:);...
                        -new_bvals.*new_bvecs(2,:).*new_bvecs(2,:);...
                        -2*new_bvals.*new_bvecs(2,:).*new_bvecs(3,:);...
                        -new_bvals.*new_bvecs(3,:).*new_bvecs(3,:)]';

                    dwv = squeeze(dw(i,j,k,:));

                    % Ignore any zero or negative measurements.
                    ignore = find(dwv<=0);
                    dwv = dwv(setdiff(1:length(dwv),ignore));

                    % Only proceed if there are enough non-zero
                    % measurements and non-zero b=0s.
                    if(length(dwv)>=7 && length(setdiff(b0_inds, ignore))>0)

                        X = X(setdiff(1:length(X),ignore), :);

                        % Use weighted linear fits.
                        W = diag(dwv);
                        X = W*X;
                        Xi = pinv(X);
                        dt(i,j,k,2:end) = Xi*W*log(dwv);
                        dt(i,j,k,1) = 0;

                        % check whether the Tensor is SPD, and 
                        % mark by 1 in first entry of dt if not                 
                        T = MakeDT_Matrix( dt(i,j,k,3),...
                                           dt(i,j,k,4),...
                                           dt(i,j,k,5),...
                                           dt(i,j,k,6),...
                                           dt(i,j,k,7),...
                                           dt(i,j,k,8) );
                        eigenvecs = eigs(T);                    
                        if ( min(eigenvecs) < 1000*eps )

                            dt(i,j,k,1) = 1 ;
                            dt_mask_invalid_tensors(i,j,k) = 1 ;

                        end

                    else
                        mask(i,j,k)=0;
                    end
                end

            end
        end
    end

    % Save each DT component as a separate nifti file.
%     dt = flipdim(dt, 1);
    for i=1:8
        write_hdr_nii(dt(:,:,:,i),...
                      [output_folder dt_pref num2str(i)], hdr);
    end


    % ------- Step 2: Compute the DTI of downsampled (low-resolution) DWI ----
    % Create low-res DTI after averaging the DW data over various nxnxn blocks 
     fprintf('Downsampling by factor of %i in each dimension (%s)\n', ...
                ds_rate, data_folders{fi});
   
    dtl = zeros(XSIZE,YSIZE,ZSIZE,8);
    dtl(:,:,:,1) = -1;
    dtl_mask_invalid_tensors = zeros(XSIZE,YSIZE,ZSIZE);

    for k=1:(ZSIZE-ds_rate+1)
        display(sprintf('Slice %i of %i (%s) ds=%i.', k, ZSIZE, data_folders{fi}, ds_rate));
        for j=1:(YSIZE-ds_rate+1)
            for i=1:(XSIZE-ds_rate+1)
                % We consider here the block of dsxdsxds with i,j,k at the corner.
                % Only process if whole block is within the mask.
                bmask = mask(i:(i+ds_rate-1),j:(j+ds_rate-1),k:(k+ds_rate-1));
                if(min(bmask(:)))

                    % Average the gradient correction over the block
                    L = reshape(squeeze(mean(mean(mean(g(i:(i+ds_rate-1),...
                                                         j:(j+ds_rate-1),...
                                                         k:(k+ds_rate-1),:)...
                                                         )))),3,3);
                    I = eye(3);

                    % correct bvecs and calculate their norm
                    v = (I+L)*bvecs;
                    n = sqrt(sum(v.^2));

                    % rotate bvecs to voxel coord system 
                    % and correct both bvecs and bvals
                    new_bvecs = v./repmat(n,3,1);
                    new_bvals = n.^2.*bvals;

                    % Construct the design matrix for linear estimation of the DT.
                    X = [ones(1, length(new_bvals));...
                        -new_bvals.*new_bvecs(1,:).*new_bvecs(1,:);...
                        -2*new_bvals.*new_bvecs(1,:).*new_bvecs(2,:);...
                        -2*new_bvals.*new_bvecs(1,:).*new_bvecs(3,:);...
                        -new_bvals.*new_bvecs(2,:).*new_bvecs(2,:);...
                        -2*new_bvals.*new_bvecs(2,:).*new_bvecs(3,:);...
                        -new_bvals.*new_bvecs(3,:).*new_bvecs(3,:)]';

                    dwv = squeeze(mean(mean(mean(dw(i:(i+ds_rate-1),...
                                                    j:(j+ds_rate-1),...
                                                    k:(k+ds_rate-1),:)...
                                                    ))));

                    % Ignore any zero or negative measurements.
                    ignore = find(dwv<=0);
                    dwv = dwv(setdiff(1:length(dwv),ignore));
                    X = X(setdiff(1:length(X),ignore), :);

                    % Use weighted linear fits.
                    W = diag(dwv);
                    X = W*X;
                    Xi = pinv(X);
                    dtl(i,j,k,2:end) = Xi*W*log(dwv);
                    dtl(i,j,k,1) = 0;

                    % check whether the Tensor is SPD, and 
                    % mark by 1 in first entry of dt if not                 
                    T = MakeDT_Matrix( dtl(i,j,k,3),...
                                       dtl(i,j,k,4),...
                                       dtl(i,j,k,5),...
                                       dtl(i,j,k,6),...
                                       dtl(i,j,k,7),...
                                       dtl(i,j,k,8) );
                    eigenvecs = eigs(T);                    
                    if ( min(eigenvecs) < 1000*eps )

                        dtl(i,j,k,1) = 1 ;
                        dtl_mask_invalid_tensors(i,j,k) = 1 ;

                    end

                end
            end
        end
    end

    % Save each DT component as a separate nifti file.
%     dtl = flipdim(dtl, 1);
    for i=1:8
        dt_lowres_name = sprintf('%slowres_%i_', dt_pref, ds_rate);
        write_hdr_nii(dtl(:,:,:,i) , ...
                     [output_folder dt_lowres_name num2str(i) ], ...
                      hdr);
    end
end

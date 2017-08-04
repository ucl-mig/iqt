function compute_dti(input_dir, output_dir, data_folders, sub_path, ...
                     dw_file, bvals_file, bvecs_file, mask_file, ...
                     dt_pref, b_max, b0_thresh)
% COMPUTE_DTI
%   Computes DTIs on the DWIs.
%
%   Args: 
%       INPUT_DIR: Input root folder (eg HCP root)
%       OUTPUT_DIR: Output root folder (could be HCP root)
%       DATA_FOLDERS: Subjects to process (eg HCP subjects)
%       SUB_PATH: Input/Output data access is done as:
%                  input_dir/data_folder/sub_path/file_name
%                  output_dir/data_folder/sub_path/file_name
%       DW_FILE: DWI data file
%       BVALS_FILE: b-values file
%       BVECS_FILE: b-vectors file
%       MASK_FILE: mask file
%       DT_PREF: DTI output prefix
%       B_MAX: max b-value to use for DTI estimation (default 1500)
%       B0_THRESH: b-value smaller than this is considered b0 (default 100)
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

if ~exist('b_max', 'var')
    b_max = 1500;
end
if ~exist('b0_thresh', 'var')
    b0_thresh = 100;
end

parfor fi = 1:length(data_folders)
    if(~exist([output_dir data_folders{fi}  '/' sub_path], 'dir'))
        mkdir(output_dir, [data_folders{fi}  '/' sub_path]);
    end
    
    fprintf('DTI computation subject: %s\n', data_folders{fi});
    input_folder = [input_dir data_folders{fi} '/' sub_path];
    output_folder = [output_dir data_folders{fi} '/' sub_path];
    
    % Read in the DWI data
    fprintf('Loading DWI: %s (%s)\n', dw_file, data_folders{fi});
    nii = load_nii( [input_folder dw_file] );
    dw = double(nii.img);
    [XSIZE,YSIZE,ZSIZE,~] = size(dw);
    hdr = nii.hdr;
    
    % Read bvals and bvecs text files
    fprintf('Loading bvals/bvecs: %s, %s (%s)\n', bvals_file, bvecs_file, data_folders{fi});
    bvecs = load([input_folder bvecs_file]); % should be 3xN
    bvals = load([input_folder bvals_file]); % should be 1xN

    % Read the brain mask
    fprintf('Loading mask: %s (%s)\n', mask_file, data_folders{fi});
    nii = load_nii( [input_folder mask_file] );
    mask = nii.img;
    nii = []; % clear to run within parfor.

    % Here we want to fit only to the inner, b=1000, shell as well as b=0
    % images.
    dw_inds = find(bvals<b_max);
    bvals = bvals(dw_inds);
    bvecs = bvecs(:,dw_inds);
    dw = dw(:,:,:,dw_inds);

    % Identify b=0 indices
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
    for i=1:8
        write_hdr_nii(dt(:,:,:,i),...
                      [output_folder dt_pref num2str(i)], hdr);
    end

end

%% Compute DTIs on the original DW images and its downsampled version.
% As a result, we obtain a pair of high-res and low-res DTIs. 

function compute_dti(input_dir, output_dir, data_folders, settings)

% Fetch the parameters:
dt_name = settings.dt_name;
ds = settings.upsample_rate; % downsampling rate

for fi = 1:length(data_folders)
    if(~exist([output_dir '/'  data_folders{fi} ]))
        mkdir([output_dir '/'  data_folders{fi} ]);
    end
    
    %% -------- Step 1: Compute DTI for the original high-res DWI ---------
    input_folder = [input_dir '/' data_folders{fi}];
    output_folder = [output_dir '/' data_folders{fi}];

    %----------------------------------------------------------------------
    % Read gradient nonlinearity file and qform (HCP only).
    filename     = 'grad_dev.nii';  % gradient nonlinearity file name
    tmp = load_nii( [input_folder '/' filename] );
    elsp = tmp.hdr.dime.pixdim(2:4); % element spacing
    g = flipdim(tmp.img,1);
    tmp = []; % clear to run within parfor.

    %----------------------------------------------------------------------
    % Read bvals and bvecs text files
    bvecs = load([input_folder '/bvecs']); % should be 3xN
    bvals = load([input_folder '/bvals']); % should be 1xN

    %----------------------------------------------------------------------
    % Read in the image data
    filename     = 'data.nii';  % gradient nonlinearity file name

    disp(['Loading NII: ' data_folders{fi}]);
    tic
    tmp = load_nii( [input_folder '/' filename] );
    toc

    disp('Flipdim');
    tic
    dw = flipdim(tmp.img,1);
    toc

    %clear tmp
    tmp = [];
    [XSIZE,YSIZE,ZSIZE,COMP] = size(dw);

    %----------------------------------------------------------------------
    % Read the brain mask
    filename     = 'nodif_brain_mask.nii';  % gradient nonlinearity file
    tmp = load_nii( [input_folder '/' filename] );
    mask = flipdim(tmp.img,1);
    %clear tmp
    tmp = [];

    % Here we want to fit only to the inner, b=1000, shell as well as b=0
    % images.
    dw_inds = find(bvals<1500);
    bvals = bvals(dw_inds);
    bvecs = bvecs(:,dw_inds);
    dw = dw(:,:,:,dw_inds);

    % Identify b=0 indices
    b0_thresh = 100;
    b0_inds = find(bvals<b0_thresh);

    %----------------------------------------------------------------------
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
        display(sprintf('Slice %i of %i (%i).', k, ZSIZE, fi));
        for j=1:YSIZE
            for i=1:XSIZE
                if(mask(i,j,k))    
                    new_bvecs = bvecs;
                    new_bvals = bvals;

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
        write_std_nii(dt(:,:,:,i),...
                      [output_folder '/' dt_name num2str(i) ], elsp );
    end

    %% ------- Step 2: Compute the DTI of downsampled (low-resolution) DWI ----
    % Create low-res DTI after averaging the DW data over various nxnxn blocks 
    
    ds = settings.upsample_rate; % up-sample rate.
    dtl = zeros(XSIZE,YSIZE,ZSIZE,8);
    dtl(:,:,:,1) = -1;
    dtl_mask_invalid_tensors = zeros(XSIZE,YSIZE,ZSIZE);

    display(sprintf('Downsampling by factor of %i in each dimension.', ds));
    for k=1:(ZSIZE-ds+1)
        display(sprintf('Slice %i of %i (%i) ds=%i.', k, ZSIZE, fi, ds));
        for j=1:(YSIZE-ds+1)
            for i=1:(XSIZE-ds+1)
                % We consider here the block of dsxdsxds with i,j,k at the corner.
                % Only process if whole block is within the mask.
                bmask = mask(i:(i+ds-1),j:(j+ds-1),k:(k+ds-1));
                if(min(bmask(:)))

                    % Average the gradient correction over the block
                    L = reshape(squeeze(mean(mean(mean(g(i:(i+ds-1),...
                                                         j:(j+ds-1),...
                                                         k:(k+ds-1),:)...
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

                    dwv = squeeze(mean(mean(mean(dw(i:(i+ds-1),...
                                                    j:(j+ds-1),...
                                                    k:(k+ds-1),:)...
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
    for i=1:8
        dt_lowres_name = sprintf('%slowres_%i_', dt_name, ds);
        write_std_nii(dtl(:,:,:,i) , ...
                     [output_folder '/' dt_lowres_name num2str(i) ], ...
                      elsp);
    end
end

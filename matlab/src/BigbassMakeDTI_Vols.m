%% Load in an HCP diffusion data set

addpath(genpath('~/Sources/UCL/iqt'));

DataPath = '~/Data/HCP';
OutputPath = '~/Data/HCP';

% Original training and test set.
%data_folders = {'101915','106016','120111','122317','130316','148335','153025','159340', '162733', '163129', '178950', '188347', '189450', '199655', '211720', '280739'};
% Additional 8 random training set elements.
%data_folders = {'351938',  '390645',  '545345',  '586460',  '705341',  '749361',  '765056',  '951457'};
% 16 22-27 white non-hispanic right-handed females.
%data_folders = {'153025',  '205725',  '103414',  '351938',  '727654',  '748258',  '151223',  '162733', '167743', '586460', '100307', '111009', '130316', '255639', '113215', '872158'};
% 8 31-36 non-white males.
%data_folders = {'182840',  '205119',  '209733',  '160123',  '528446',  '581349',  '165840',  '211417'};
% Diverse training and testing set.
%data_folders = {'992774', '904044', '165840', '654754', '125525', '889579', '205119', '713239', '570243', '133928', '899885', '448347', '117324', '214423', '857263', '153025'};
% Lots of other ones for testing variance.
% data_folders = {'106319', '117122', '133827', '140824', '158540', '196750', '205826', '366446', '685058', '734045', '826353', '887373',...
% '100408', '110411', '126325', '142828', '151526', '159239', '169343', '197550', '208226', '217126', '371843', '530635', '598568', '688569', '856766', '959574',...
% '118730', '127933', '134324', '143325', '151627', '172332', '190031', '198451', '217429', '284646', '541943', '627549', '702133', '894673', '978578',...
% '102816', '111312', '118932', '128632', '135932', '144226', '175439', '191437', '199150', '209935', '221319', '293748', '397760', '638049', '704238', '751348', '859671', '896879', '984472',...
% '111514', '119833', '129028', '136833', '148032', '153429', '161731', '176542', '192439', '210617', '224022', '298051', '414229', '547046', '645551', '753251', '861456', '899885',...
% '103515', '111716', '130013', '137128', '154431', '162329', '177746', '192540', '200614', '239944', '304020', '429040', '559053', '756055', '865363', '901139',...
% '103818', '112819', '120212', '130316', '138231', '149337', '156233', '193239', '201111', '245333', '307127', '561242', '665254', '715647', '761957', '871964',...
% '105115', '130922', '138534', '149539', '156637', '182739', '194140', '201414', '212318', '246133', '329440', '485757', '672756', '917255',...
% '105216', '113619', '123117', '131924', '139637', '150423', '157336', '163432', '195647', '214019', '249947', '497865', '579665', '677968', '729557', '788876', '877168', '932554',...
% '115320', '124422', '133625', '140420', '150524', '158035', '185139', '196144', '214221', '250427', '355239', '499566', '680957', '732243', '792564', '885975', '937160'};
data_folders = {'100307'};


b0_thresh = 100


for fi = 1:length(data_folders)
  if(~exist([OutputPath '/'  data_folders{fi} ]))
    mkdir([OutputPath '/'  data_folders{fi} ]);
    mkdir([OutputPath '/'  data_folders{fi} '/T1w' ]);
    mkdir([OutputPath '/'  data_folders{fi} '/T1w/Diffusion/' ]);
  end
    
    data_folder = [ DataPath '/' data_folders{fi} '/T1w/Diffusion' ];

    output_dir = [OutputPath '/'  data_folders{fi} '/T1w/Diffusion' ];
    
    
    %--------------------------------------------------------------------------
    % Read gradient nonlinearity file and qform
    filename     = 'grad_dev.nii';  % gradient nonlinearity file name
    tmp = load_nii( [data_folder '/' filename] );
    elsp = tmp.hdr.dime.pixdim(2:4); % element spacing
    g = flipdim(tmp.img,1);
    %clear tmp
    % Need to do this instead to run within parfor.
    tmp = [];


    %--------------------------------------------------------------------------
    % Read bvals and bvecs text files
    bvecs = load([data_folder '/bvecs']); % should be 3xN
    bvals = load([data_folder '/bvals']); % should be 1xN


    %--------------------------------------------------------------------------
    % Read in the image data
    filename     = 'data.nii';  % gradient nonlinearity file name

    disp(['Loading NII: ' data_folders{fi}]);
    tic
    tmp = load_nii( [data_folder '/' filename] );
    toc

    disp('Flipdim');
    tic
    dw = flipdim(tmp.img,1);
    toc

    %clear tmp
    tmp = [];

    [XSIZE,YSIZE,ZSIZE,COMP] = size(dw);


    %--------------------------------------------------------------------------
    % Read the brain mask
    filename     = 'nodif_brain_mask.nii';  % gradient nonlinearity file name
    tmp = load_nii( [data_folder '/' filename] );
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
    b0_inds = find(bvals<b0_thresh);
    
    %--------------------------------------------------------------------------
    % do the high-resolution reconstruction
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
                    
                    % create matrices
%                     L = reshape(squeeze(g(i,j,k,:)),3,3);
%                     I = eye(3);
%     
%                     % correct bvecs and calculate their norm
%                     v = (I+L)*bvecs;
%                     n = sqrt(sum(v.^2));
%     
%                     % rotate bvecs to voxel coord system 
%                     % and correct both bvecs and bvals
%                     new_bvecs = v./repmat(n,3,1);
%                     new_bvals = n.^2.*bvals;
    
                    new_bvecs = bvecs;
                    new_bvals = bvals;
                    
                    % Construct the design matrix for linear estimation of the DT.
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

                        % check whether the Tensor is SPD, and mark by 1 in first entry of dt if not                 
                        T = MakeDT_Matrix( dt(i,j,k,3),dt(i,j,k,4),dt(i,j,k,5),dt(i,j,k,6),dt(i,j,k,7),dt(i,j,k,8) );
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
    
    for i=1:8
        write_std_nii( dt(:,:,:,i) , [output_dir '/dt_b1000_' num2str(i) ] , elsp );
    end
    %write_std_nii( dt_mask_invalid_tensors , 'dt_mask_invalid_tensors_b1000' , elsp );
    
    %--------------------------------------------------------------------------
    % Create low resolution diffusion tensor image after averaging the dw data over various nxnxn blocks 
    for ds = 2:4

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
                        L = reshape(squeeze(mean(mean(mean(g(i:(i+ds-1),j:(j+ds-1),k:(k+ds-1),:))))),3,3);
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

                        dwv = squeeze(mean(mean(mean(dw(i:(i+ds-1),j:(j+ds-1),k:(k+ds-1),:)))));

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

                        % check whether the Tensor is SPD, and mark by 1 in first entry of dt if not                 
                        T = MakeDT_Matrix( dtl(i,j,k,3),dtl(i,j,k,4),dtl(i,j,k,5),dtl(i,j,k,6),dtl(i,j,k,7),dtl(i,j,k,8) );
                        eigenvecs = eigs(T);                    
                        if ( min(eigenvecs) < 1000*eps )

                            dtl(i,j,k,1) = 1 ;
                            dtl_mask_invalid_tensors(i,j,k) = 1 ;

                        end

                    end
                end
            end
        end

        for i=1:8
            write_std_nii( dtl(:,:,:,i) , [ output_dir '/dt_b1000_lowres_' num2str(ds) '_' num2str(i) ] , elsp );
        end
        %write_std_nii( dtl_mask_invalid_tensors , ['dtl_mask_invalid_tensors_b1000' num2str(ds)] , elsp );
    end
    
   
    cd('../../..');
        
end



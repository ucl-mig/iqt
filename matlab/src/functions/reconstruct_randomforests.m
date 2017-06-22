function reconstruct_randomforests(input_dir, output_dir, trees_dir, ...
                                   data_folders, sub_path, dt_pref, ...
                                   ds_rate, no_rnds, sample_rate, ...
                                   input_radius, fv, edge_recon, ...
                                   flip_dim)
% RECONSTRUCT_RANDOMFORESTS
%   Computes IQT super-resolution DTI from the input DTI (low-resolution).
%
%   Args: 
%       INPUT_DIR: Input root folder of low-res DTI
%       OUTPUT_DIR: Output root folder
%       TREES_DIR: Dir containing the RF trees
%       DATA_FOLDERS: Subjects to process (eg HCP subjects)
%       SUB_PATH: Input/Output data access is done as:
%                  input_dir/data_folder/sub_path/file_name
%                  output_dir/data_folder/sub_path/file_name
%       DT_PREF: DTI output prefix
%       DS_RATE: Super-resolution factor
%       NO_RNDS: number of trees in the RF
%       SAMPLE_RATE: data sampling rate for random sub-sampling
%       INPUT_RADIUS: the input is a cubic patch of size (2*INPUT_RADIUS+1)^3
%       FV: feature version used for computing features for tree training
%       EDGE_RECON: flag that indicates whether to reconstructs edge or not
%                   (edge reconstruction is slow!)
%       FLIP_DIM (optional): The internal IQT code flips all images along 
%                            axis-1, due to historic reasons. But this 
%                            should not be the case for user data. 
%                            Normally, this flag should be off (default 0).
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
check_path(trees_dir);
check_path(sub_path);

rescale_factor = 1.0;
scale_const = 1E-3;

ds=ds_rate;
n=input_radius;
m=ds;

if ~exist('flip_dim', 'var')
    flip_dim = 0;
end

% Trees:   
tail_name = sprintf('DS%02i_%ix%ix%i_%ix%ix%i_Sub%03i', ds, 2*n+1,2*n+1,2*n+1, m, m, m, sample_rate);
tree_name = ['RegTreeValV' int2str(fv) '_' tail_name '_'];

trees=cell(no_rnds,1); 
for k= 1:no_rnds
    tmp=load(sprintf([trees_dir tree_name '%04i.mat'], k));
    trees{k} = tmp.tree;
end


% Perform reconstruction: 
parfor dataid = 1:length(data_folders)
    fprintf(['\nReconstructing: ' data_folders{dataid} '\n']);
    
    output_folder = [output_dir data_folders{dataid} '/'];
    if(~exist(output_folder, 'dir'))
        mkdir(output_folder);
    end
   
    
    % Load in the diffusion tensor images (low-res)
    [dt_lr, hdr] = ReadDT_Volume([ input_dir data_folders{dataid} '/' sub_path dt_pref]);
    x_h=[]; y_h=[]; z_h=[];
    %is this typically HCP data?
    if flip_dim>0
        % Load in hires diffusion tensor images
        dth = dt_lr;
        [x_h, y_h, z_h, ~] = size(dth);
        dth = [];
        % low res HCP
        dt_lr = ReadDT_Volume([ input_dir data_folders{dataid}  '/' sub_path dt_pref sprintf('lowres_%i_', ds)]);
        dto = dt_lr(1:ds:end,1:ds:end,1:ds:end,:); % iput low-res dti.
    else
        dto = dt_lr;
    end

    % Reconstruct:
    if edge_recon % Reconstruct on the boundary.
        % compute the features:
        FeatureMapFile = [input_dir data_folders{dataid} '/' sub_path 'FeaturesMapEdge_V' int2str(fv) '_' tail_name '.mat'];
        if(exist(FeatureMapFile, 'file'))
            display(['Loading feature map: ' FeatureMapFile])
            data = load(FeatureMapFile);
            dto_features = data.dto_features;
            T = data.T;
            comipatchmean = data.comipatchmean;
            data = [];
        else
            fprintf('Computing feature map... (%s)\n', data_folders{dataid});
            patchlib_name = [patchlibs_dir '/PatchLibs_' tail_name '_0001.mat'];
            data = load(patchlib_name);
            comipatchlib = data.comipatchlib;
            data = [];
            comipatchmean = mean(comipatchlib);
            ncomipatchlib = comipatchlib - repmat(comipatchmean, [length(comipatchlib), 1]);
            T = ncomipatchlib'*ncomipatchlib/length(ncomipatchlib);
            dto_features = PatchFeatureMapEdge(dto, n, T, comipatchmean, fv);
            myS = [];
            myS.dto_features = dto_features; myS.T = T; myS.comipatchmean = comipatchmean;
            parsave_struct(FeatureMapFile, myS, '-v7.3');
            myS = [];
        end

        % reconstruct and save the estimated DTI and precision map:
       output_subdir  = sprintf(['RF_Edge_V' int2str(fv) '_NoTree%02i_' tail_name '/'], no_rnds);
       [img_recon, img_confid] = ForestSuperResEdge(dto, trees, n, m, dto_features, ds, T, comipatchmean, rescale_factor, scale_const, 'weighted_average');
    else % perform super-resolution on the interior.
        % compute the features:
        FeatureMapFile = [input_dir data_folders{dataid} '/' sub_path 'Features_V' int2str(fv) '_' tail_name '.mat'];
        if(exist(FeatureMapFile, 'file'))
            display(['Loading feature map: ' FeatureMapFile])
            data = load(FeatureMapFile);
            featuresMap = data.featuresMap;
            data = [];
        else
            fprintf('Computing feature map... (%s)\n', data_folders{dataid});
            featuresMap = PatchFeatureMap(dto, n, m, fv);
            parsave(FeatureMapFile, featuresMap);
        end

        % reconstruct and save the estimated DTI and precision map:
       output_subdir  = sprintf(['RF_V' int2str(fv) '_NoTree%02i_' tail_name '/'], no_rnds);
       ValidVar = 1;
       [img_recon, img_confid, img_var] = ForestSuperRes(dto, trees, n, m, featuresMap, ds, rescale_factor, ValidVar);    
    end

    %is this typically HCP data?
    if flip_dim>0
        % crop so the output size is the same as the ground truth.
        img_recon = img_recon(1:x_h,1:y_h,1:z_h,:);
    else
        hdr.dime.pixdim(2:4) = hdr.dime.pixdim(2:4) / ds;
        img_recon = flipdim(img_recon, 1);
    end
    
    % Save all reconstructed files:
    recon_dir = [output_folder sub_path output_subdir];
    fprintf('Saving outputs to: %s (%s)\n', recon_dir, data_folders{dataid});
    if(~exist(recon_dir, 'dir'))
        mkdir(recon_dir);
    end
    %parsave([output_folder output_subdir 'settings.mat'],'settings') ;
    
    % Store the fitted DTI as mat files named dt_recon_i, i = 1,...,8
    for i=1:8 
        write_hdr_nii( img_recon(:,:,:,i), [recon_dir dt_pref 'recon_' num2str(i)], hdr );
    end
    
    % Store the precision (1/variance) over the estimated DTI as nifti
    % files.
    write_hdr_nii(img_confid, [recon_dir dt_pref 'confid'], hdr);
    
end





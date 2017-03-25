function reconstruct_randomforests(data_folders, settings)

%% Configurations:
% Paths:
input_dir = settings.input_dir;
output_dir = settings.output_dir;
trees_dir = settings.trees_dir;
trees_list = settings.trees_list;
patchlibs_dir = settings.patchlibs_dir;
dt_name = settings.dt_name;
rescale_factor = 1.0;
scale_const = 1E-3;

% Data:
subsample_rate=settings.subsample_rate;
fv=settings.feature_version;

% Problem:
ds=settings.upsample_rate;
n=settings.input_radius;
m=ds;
edge_recon = settings.edge; % true if you want to reconstruct on the edge

% Trees:   
tail_name = sprintf('DS%02i_%ix%ix%i_%ix%ix%i_Sub%03i', ds, 2*n+1,2*n+1,2*n+1, m, m, m, subsample_rate)
tree_name = ['RegTreeValV' int2str(fv) '_' tail_name '_'];

trees={}; 
for k= 1:length(trees_list)
    tmp=load(sprintf([trees_dir '/' tree_name '%04i.mat'], trees_list(k)));
    trees{k} = tmp.tree;
end


%% Perform reconstruction: 
for dataid = 1:length(data_folders)
    display(sprintf(['\nReconstructing: ' data_folders{dataid} '\n']))
    
    output_folder = [output_dir '/' data_folders{dataid}];
    if(~exist(output_folder))
        mkdir(output_folder);
    end
   
    
    % Load in the original and low-res diffusion tensor images:
    % original res
    file_orig = [input_dir '/' data_folders{dataid} '/' dt_name];
    dth = ReadDT_Volume(file_orig);
    [x_h, y_h, z_h, junk] = size(dth);
    clear dth
    
    % low-res:
    dt_lr = ReadDT_Volume([ input_dir '/' data_folders{dataid} '/' dt_name sprintf('lowres_%i_', ds)]);
    dto = dt_lr(1:ds:end,1:ds:end,1:ds:end,:); % iput low-res dti.

    % Reconstruct:
    if edge_recon % Reconstruct on the boundary.
        % compute the features:
        FeatureMapFile = [input_dir '/' data_folders{dataid} '/FeaturesMapEdge_V' int2str(fv) '_' tail_name '.mat'];
        if(exist(FeatureMapFile))
            load(FeatureMapFile)
        else
            display('Computing feature map...');
            patchlib_name = [patchlibs_dir '/PatchLibs_' tail_name '_0001.mat'];
            display(['Loading the patch library: ' patchlib_name])
            load(patchlib_name)
            comipatchmean = mean(comipatchlib);
            ncomipatchlib = comipatchlib - repmat(comipatchmean, [length(comipatchlib), 1]);
            T = ncomipatchlib'*ncomipatchlib/length(ncomipatchlib);
            dto_features = PatchFeatureMapEdge(dto, n, T, comipatchmean, fv);
            save(FeatureMapFile, 'dto_features', 'T', 'comipatchmean');
        end

        % reconstruct and save the estimated DTI and precision map:
       output_subdir  = sprintf(['RF_Edge_V' int2str(fv) '_NoTree%02i_' tail_name], length(trees_list));
       [img_recon, img_confid] = ForestSuperResEdge(dto, trees, n, m, dto_features, ds, T, comipatchmean, rescale_factor, scale_const, 'weighted_average');
    else % perform super-resolution on the interior.
        % compute the features:
        FeatureMapFile = [input_dir '/' data_folders{dataid} '/Features_V' int2str(fv) '_' tail_name '.mat'];
        if(exist(FeatureMapFile))
            display(['Loading feature map: ' FeatureMapFile])
            load(FeatureMapFile)
        else
            display('Computing feature map...');
            featuresMap = PatchFeatureMap(dto, n, m, fv);
            save(FeatureMapFile, 'featuresMap');
        end

        % reconstruct and save the estimated DTI and precision map:
       output_subdir  = sprintf(['RF_V' int2str(fv) '_NoTree%02i_' tail_name], length(trees_list));
       ValidVar = 1;
       [img_recon, img_confid, img_var] = ForestSuperRes(dto, trees, n, m, featuresMap, ds, rescale_factor, ValidVar);    
    end

    % crop so the output size is the same as the ground truth.
    img_recon = img_recon(1:x_h,1:y_h,1:z_h,:);
    
    % Save all reconstructed files:
    display(['Saving outputs to: ' output_folder '/' output_subdir])
    if(~exist([output_folder '/' output_subdir]))
        mkdir([output_folder '/' output_subdir]);
    end
    %save([output_folder '/' output_subdir '/dt_recon.mat'], 'img_recon');
    %save([output_folder '/' output_subdir '/dt_confid.mat'],'img_confid') ;
    save([output_folder '/' output_subdir '/settings.mat'],'settings') ;
    
    % Store the fitted DTI as mat files named dt_b1000_i, i = 1,...,8
    for i=1:8 
        tic
        write_std_nii( img_recon(:,:,:,i), [output_folder '/' output_subdir '/dt_recon_' num2str(i) ]);
        toc
    end
    
    % Store the precision (1/variance) over the estimated DTI as nifti
    % files.
     for i=1:8 
        tic
        write_std_nii( img_confid(:,:,:,i), [output_folder '/' output_subdir '/dt_confid_' num2str(i) ]);
        toc
    end
   
end





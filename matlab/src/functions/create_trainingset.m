function create_trainingset(input_dir, train_dir, data_folders, ...
                            sub_path, sample_rate, no_rnds, ds_rate, input_radius)
% CREATE_TRAININGSET samples the list of full patch-pairs for the training
%   subjects and sub-samples randomly to build multiple training patch
%   libraries that are smaller in size to fit within the RAM constraints.
%
%   Args:
%       INPUT_DIR: Input root folder (eg HCP root)
%       TRAIN_DIR: This directory will contain all training outputs
%       DATA_FOLDERS: Subjects to process (eg HCP subjects)
%       SUB_PATH: Input/Output data access is done as:
%                  input_dir/data_folder/sub_path/file_name
%                  output_dir/data_folder/sub_path/file_name
%       SAMPLE_RATE: data sampling rate for random sub-sampling
%       NO_RNDS: number of random sub-samples to create these many training
%                datasets
%       DS_RATE: Super-resolution factor
%       INPUT_RADIUS: the input is a cubic patch of size (2*INPUT_RADIUS+1)^3
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


SR = sample_rate;
rnds = no_rnds;
n = input_radius; % the radius of the low-res patch. 
ds = ds_rate;
m = ds;

if(~exist(train_dir, 'dir'))
    mkdir(train_dir);
end

% Loop over randomizations
for si = 1:rnds
    for i=1:length(data_folders)
        load([input_dir '/' data_folders{i} '/' ...
              sprintf('ipatchlibDS%02i_N%02i', ds, n)]);
        load([input_dir '/' data_folders{i} '/' ...
              sprintf('opatchlibDS%02i_N%02i_M%02i', ds, n, m)]);
          
%         if(n>1)
%             % Need to subsample the opatchlib to match the ipatchlib
%             load([intput_dir '/' data_folders{i} sprintf('indices_to_N1_indicesDS%02i_N%02i', ds, n)]);
%             opatchlib = opatchlib(indices_to_N1_indices,:);
%         end

        keep = find(rand(1,length(ipatchlib))<(1/SR));
        if(i==1)
            comipatchlib = ipatchlib(keep,:);
            comopatchlib = opatchlib(keep,:);
        else
            comipatchlib = [comipatchlib; ipatchlib(keep,:)];
            comopatchlib = [comopatchlib; opatchlib(keep,:)];
        end
        
        % Keep track of which voxels are included.
        patchlibindices{i} = keep;
    end

    % Save some space...
    comipatchlib = single(comipatchlib);
    comopatchlib = single(comopatchlib);

    filename = sprintf('PatchLibs_DS%02i_%ix%ix%i_%ix%ix%i_Sub%03i_%04i.mat', ds, 2*n+1,2*n+1,2*n+1, m, m, m, SR, si);
    fprintf('Saving the patch-library as: %s\n', filename);
    save([train_dir filename], 'comipatchlib', 'comopatchlib', 'patchlibindices', '-v7.3');
end
fprintf('The training set is available at: %s\n', train_dir);



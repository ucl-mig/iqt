%% Reconstruction:

%Header of DT images.
settings.dt_name = 'dt_b1000_';
%Radius of the low-res input patch.
settings.input_radius = 2;
%Upsampling rate
settings.upsample_rate = 2; 
%Subsampling rate.
settings.subsample_rate = 32; 
%Feature set used in Neuroimage paper. See PatchFeatureList.m for details.
settings.feature_version = 6; 
% Set it to 1 to perform boundary completion. 
settings.edge = 0;

% Paths:
settings.input_dir = '/dir/where/input/DTIs/are/stored';
settings.output_dir = '/dir/where/output/DTIs/are/saved';
settings.trees_dir = '/dir/where/input/trees/are/stored';
settings.trees_list = [1:8]; % indixes of trees included in the forest.
settings.patchlibs_dir = '/dir/where/training/patch-libs/are/stored';
data_folders = {'subject_1','subject_2','subject_3'}; % list of sub-dirs.

% Perform super-resolution: 
reconstruct_randomforests(data_folders, settings)

% Visualise the MD/FA/CFA: (?) save as a .png?
visualise_results(data_folders, settings)






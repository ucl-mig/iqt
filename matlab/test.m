% TEST  A script that uses pretrained trees that are provided with
%   the codebase to reconstruct typical HCP low resolution DTI.
%
%   This is a script, you have to edit the 'settings'.
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

%% Settings
addpath(genpath('.'));

%Header of DT images.
settings.dt_name = 'dt_b1000_';
%Radius of the low-res input patch i.e. the input is a cubic patch of size (2*settings.input_radius+1)^3
settings.input_radius = 2;
%Upsampling rate
settings.upsample_rate = 2; 
%Subsampling rate. This determines the size of training sets. From each subject, we randomly draw patches with probability 1/settings.subsample_rate 
settings.subsample_rate = 32; 
%Feature set used in Neuroimage paper. They are a combination of scale/rotationally invariant features. See PatchFeatureList.m for details.
settings.feature_version = 6; 
% Set it to 1 to perform boundary completion. By default, set to 0 and
% performs reconstruction only on the interior region of the brain.
settings.edge = 0;

% Paths:
settings.input_dir = '/dir/where/input/DTIs/are/stored';
settings.output_dir = '/dir/where/output/DTIs/are/saved';
settings.trees_dir = '/dir/where/input/trees/are/stored';
settings.trees_list = [1:8]; % indices of trees included in the forest.
settings.patchlibs_dir = '/dir/where/training/patch-libs/are/stored';
data_folders = {'subject_1','subject_2','subject_3'}; % list of sub-dirs.


%% Reconstruction
% Perform super-resolution: 
reconstruct_randomforests(data_folders, settings)


%% Visualisation
% Visualise the MD/FA/CFA ans save the figure as a FIG file.
visualise_results(data_folders, settings)






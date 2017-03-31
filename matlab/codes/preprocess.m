%% Perform preprocessing and generate training data:

% Set the paths:
dwi_dir = '/dir/where/DWIs/are/stored'; 
dti_dir = '/dir/where/DTIs/are/saved';
traindata_dir = '/dir/where/training/sets/are/saved';
data_folders = {'subject_1','subject_2','subject_3'}; % list of sub-dirs.

settings.dt_name = 'dt_b1000_';
settings.input_radius = 2; % the radius of the low-res input patch i.e. the input is a cubic patch of size (2*settings.input_radius+1)^3
settings.upsample_rate = 2; % the upsampling rate
settings.subsample_rate = 32; % determines the size of training sets. From each subject, we randomly draw patches with probability 1/settings.subsample_rate 
settings.no_rnds = 8; % no of separate training sets to be created. 

% Step 1: model computation (e.g. DTI) from DWIs.
% the high-res and low-res DTI are computed from the DWIs and its artificially
% downsampled version. 
compute_dti(dwi_dir, dti_dir, data_folders, settings)

% Step 2: create a library of low-res and high-res patches.
% for each subject, the exhaustive list of all patch pairs are saved in a large
% matrix.
compute_patchlib(dti_dir, dti_dir, data_folders, settings)

% Step 3: create a training set. 
% from the library of each subjct, we randomly select a subset of patch
% pairs with probability 1/settings.subsample_rate.
% repeat this process settings.no_rnds number of times to create the
% specified number of separate trainign sets. 
create_trainingset(dti_dir, traindata_dir, data_folders, settings)
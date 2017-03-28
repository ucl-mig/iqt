%% Perform preprocessing and generate training data:

% Set the paths:
dwi_dir = '/dir/where/DWIs/are/stored'; 
dti_dir = '/dir/where/DTIs/are/stored';
traindata_dir = '/dir/where/training/sets/are/stored';
data_folders = {'subject_1','subject_2','subject_3'}; % list of sub-dirs.

settings.dt_name = 'dt_b1000_';
settings.input_radius = 2; % the radius of the low-res input patch.
settings.upsample_rate = 2; % the upsampling rate
settings.subsample_rate = 32; % the rate of subsampling.
settings.no_rnds = 8; % no of randomisations

% Step 1: model computation (e.g. DTI) from DWIs.
% compute the high-res and low-res DTI.
compute_dti(dwi_dir, dti_dir, data_folders, settings)

% Step 2: create a library of low-res and high-res patches.
compute_patchlib(dti_dir, dti_dir, data_folders, settings)

% Step 3: create a training set. 
create_trainingset(dti_dir, traindata_dir, data_folders, settings)
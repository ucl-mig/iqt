% IQT for DTI Super-Resolution:
% This script performs preprocessing and generates training set.

% Set the paths:
dwi_dir = '/Users/ryutarotanno/test_iqt/HCP'; 
dti_dir = data_dir;
traindata_dir = '~/tmp/iqt_codes/training_data';

params.dt_name = 'dt_b1000_';
params.input_radius = 2; % the radius of the low-res input patch.
params.upsample_rate = 2; % the upsampling rate
params.subsample_rate = 32; % the rate of subsampling.
params.no_rnds = 8; % no of randomisations


%data_folders = {'904044', '165840', '889579', '713239', ...
%                '899885', '117324', '214423', '857263'}; 
data_folders = {'117324/T1w/Diffusion'};

% Step 1: model computation (e.g. DTI) from DWIs.
% compute the high-res and low-res DTI.
compute_dti(dwi_dir, dti_dir, data_folders, params)

% Step 2: create a library of low-res and high-res patches.
compute_patchlib(dti_dir, dti_dir, data_folders, params)

% Step 3: create a training set. 
create_trainingset(dti_dir, traindata_dir, data_folders, params)
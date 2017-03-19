% Testing the pipeline locally:

%% Preprocessing:
% Set the paths:
dwi_dir = '/Users/ryutarotanno/test_iqt/HCP'; 
dti_dir = '~/tmp/iqt_codes';
traindata_dir = '~/tmp/iqt_codes/training_data';

settings.dt_name = 'dt_b1000_';
settings.input_radius = 2; % the radius of the low-res input patch.
settings.upsample_rate = 2; % the upsampling rate
settings.subsample_rate = 32; % the rate of subsampling.
settings.no_rnds = 8; % no of randomisations


%data_folders = {'904044', '165840', '889579', '713239', ...
%                '899885', '117324', '214423', '857263'}; 
data_folders = {'117324/T1w/Diffusion'};

% Step 1: model computation (e.g. DTI) from DWIs.
% compute the high-res and low-res DTI.
compute_dti(dwi_dir, dti_dir, data_folders, settings)

% Step 2: create a library of low-res and high-res patches.
compute_patchlib(dti_dir, dti_dir, data_folders, settings)

% Step 3: create a training set. 
create_trainingset(dti_dir, traindata_dir, data_folders, settings)


%% Training of trees:




%% Testing (reconstruction):


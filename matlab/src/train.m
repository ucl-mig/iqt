% Train trees:
params.dt_name = 'dt_b1000_';
params.input_radius = 2; % the radius of the low-res input patch.
params.upsample_rate = 2; % the upsampling rate
params.subsample_rate = 32; % the rate of subsampling.
params.no_rnds = 8; % no of randomisations

% Set the paths:
dwi_dir = '/Users/ryutarotanno/test_iqt/HCP'; 
dti_dir = '~/tmp/iqt_codes';
traindata_dir = '~/tmp/iqt_codes/training_data';
data_folders = {'117324/T1w/Diffusion'};



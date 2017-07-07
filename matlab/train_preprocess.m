% TRAIN_PREPROCESS  A script that creates training data for IQT random forest
%   training from a typical HCP dataset.
%   Typical usage order: TRAIN_PREPROCESS, TRAIN_RF, TEST_RF (optional)
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
% -- edit the following settings --

addpath(genpath('.'));

% Set paths (always end directory paths with a forward/back slash)
inp_dir = '/HCP/root/'; % dir where DWI data is stored (eg HCP data root)
out_dir = '/TrainResults/root/';  % typically root dir where results are stored
% list of training data subjects
data_folders = {'992774', '125525'}; %, '205119', '133928', '570243', '448347', '654754', '153025'}; 

train_dir = [out_dir 'TrainingData/']; % dir where training sets will be saved

% Check
if strcmp(inp_dir, '') || strcmp(out_dir, '')
    error('[IQT] Input/Output root missing, please check paths in settings.');
end

% Optional settings
sub_path = 'T1w/Diffusion/'; % internal directory structure
dw_file = 'data.nii'; % DWI file
bvals_file = 'bvals'; % b-values file
bvecs_file = 'bvecs'; % b-vectors file
mask_file = 'nodif_brain_mask.nii'; % mask file
grad_file = 'grad_dev.nii'; % gradient non-linearities (HCP only: grad_dev.nii)
                            % For non-HCP: grad_file = ''
dt_pref = 'dt_b1000_'; % DTI name prefix

upsample_rate = 2; % the super-resolution factor (m in paper)
input_radius = 2; % the radius of the low-res input patch i.e. the input is a cubic patch of size (2*input_radius+1)^3 (n in paper)
datasample_rate = 2; % determines the size of training sets. From each subject, we randomly draw patches with probability 1/datasample_rate
no_rnds = 8; % no of separate training sets to be created

% -- end of settings --

%%
open_matlabpool();

        
%% Step 1: model computation (e.g. DTI) from DWIs.
% the high-res and low-res DTI are computed from the DWIs by artificially
% downsampling. 
compute_dti_respairs(inp_dir, out_dir, data_folders, sub_path, ...
                     upsample_rate, dw_file, bvals_file, bvecs_file, ...
                     mask_file, grad_file, dt_pref);


%% Step 2: create a library of low-res and high-res patches.
% for each subject, the exhaustive list of all patch pairs are saved in a large
% matrix.
compute_patchlib(out_dir, out_dir, data_folders, sub_path, ...
                 dt_pref, upsample_rate, input_radius);


%% Step 3: create a training set. 
% from the library of each subjct, we randomly select a subset of patch
% pairs with probability 1/settings.subsample_rate.
% repeat this process settings.no_rnds number of times to create the
% specified number of separate trainign sets. 
create_trainingset(out_dir, train_dir, data_folders, sub_path, ...
                   datasample_rate, no_rnds, upsample_rate, input_radius);

               
%%
close_matlabpool();

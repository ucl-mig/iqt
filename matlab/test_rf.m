% TEST_RF  A (test) script.
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
inp_dir = '/HCP/root'; % dir where DWI data is stored (eg HCP data root)
out_dir = '/Results/root';  % typically root dir where results are stored
% list of test data subjects
data_folders = {'904044', '165840'}; %, '889579', '713239', '899885', '117324', '214423', '857263'};

trees_dir = './trees/'; % dir where RF trees are saved (default: precomputed trees dir)

% Check
if strcmp(out_dir, '')
    error('[IQT] Output root missing, please check paths in settings.');
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
feature_version = 6; % feature set used in Neuroimage paper. See PatchFeatureList.m for details

% Set it to 1 to perform boundary completion. By default, set to 0 and
% performs reconstruction only on the interior region of the brain.
% *** Beware: slow! (start without edge for quick reconstruction) ***
construct_edge = 0;

% The internal IQT code flips all images along axis-1, due to historic
% reasons. But this should not be the case for user data. Normally, this
% flag should be off.
flip_dim = 1;

% -- end of settings --

%%
open_matlabpool();


%% Step 1: Synthetic testing: model computation (e.g. DTI) from DWIs.
% the high-res and low-res DTI are computed from the DWIs by artificially
% downsampling. 
compute_dti_respairs(inp_dir, out_dir, data_folders, sub_path, ...
                     upsample_rate, dw_file, bvals_file, bvecs_file, ...
                     mask_file, grad_file, dt_pref);

                 %% Reconstruction
% Perform super-resolution: 
reconstruct_randomforests(out_dir, out_dir, trees_dir, data_folders, ...
                          sub_path, dt_pref, upsample_rate, no_rnds, ...
                          datasample_rate, input_radius, feature_version, ...
                          construct_edge, flip_dim);


%% Visualisation
% Visualise the MD/FA/CFA ans save the figure as a FIG file.
visualise_results(out_dir, out_dir, out_dir, data_folders, sub_path, ...
                  dt_pref, upsample_rate, no_rnds, datasample_rate, ...
                  input_radius, feature_version, construct_edge, flip_dim);

                      
%%
close_matlabpool();

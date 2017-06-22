% RECONSTRUCT_HIRES_DTI  A script to compute super-resolved DTI from the
%   pre-computed trees trained on HCP data (corresponding to the paper).
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

% Set paths (always end directory paths with a forward/back slash)
inp_dir = '/cs/research/vision/hcp/Auro/iqt.github_test/'; % root dir where DTI data is stored
out_dir = inp_dir;  % typically root dir where results are stored
train_dir = [out_dir 'TrainingData/']; % dir where training sets are saved
% list of test data subjects
data_folders = {'904044', '165840'}; %, '889579', '713239', '899885', '117324', '214423', '857263'};

% Check
if strcmp(inp_dir, '') || strcmp(out_dir, '')
    error('[IQT] Input/Output root missing, please check paths in settings.');
end

% Optional settings
sub_path = 'T1w/Diffusion/'; % internal directory structure
dt_pref = 'dt_b1000_'; % DTI name prefix

upsample_rate = 2; % the super-resolution factor
input_radius = 2; % the radius of the low-res input patch i.e. the input is a cubic patch of size (2*input_radius+1)^3
datasample_rate = 32; % determines the size of training sets. From each subject, we randomly draw patches with probability 1/datasample_rate
no_rnds = 8; % no of separate training sets to be created
feature_version = 6; % feature set used in Neuroimage paper. See PatchFeatureList.m for details

% Set it to 1 to perform boundary completion. By default, set to 0 and
% performs reconstruction only on the interior region of the brain.
% *** Beware: slow! (start without edge for quick reconstruction) ***
construct_edge = 0;

% The internal IQT code flips all images along axis-1, due to historic
% reasons. But this should not be the case for user data. Normally, this
% flag should be off.
flip_dim = 0;


%%
open_matlabpool();


%% Reconstruction
% Perform super-resolution: 
reconstruct_randomforests(inp_dir, out_dir, train_dir, data_folders, ...
                          sub_path, dt_pref, upsample_rate, no_rnds, ...
                          datasample_rate, input_radius, feature_version, ...
                          construct_edge);


%%
close_matlabpool();
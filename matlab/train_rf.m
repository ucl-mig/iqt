% TRAIN_RF  A script that trains random forest trees from training data that has
%   been prepared from a typical HCP dataset.
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
out_dir = '/TrainResults/root/'; % typically root dir where results are stored

train_dir = [out_dir 'TrainingData/']; % dir where training sets will be saved

% Check
if strcmp(out_dir, '')
    error('[IQT] Output root missing, please check paths in settings.');
end

% Optional settings
upsample_rate = 2; % the super-resolution factor (m in paper)
input_radius = 2; % the radius of the low-res input patch i.e. the input is a cubic patch of size (2*input_radius+1)^3 (n in paper)
datasample_rate = 2; % determines the size of training sets. From each subject, we randomly draw patches with probability 1/datasample_rate
no_rnds = 8; % no of separate training sets to be created
feature_version = 6; % feature set used in Neuroimage paper. See PatchFeatureList.m for details

% -- end of settings --

%%
open_matlabpool();


%% Train trees:
train_trees(train_dir, upsample_rate, input_radius, ...
            datasample_rate, no_rnds, feature_version);

        
%%
close_matlabpool();

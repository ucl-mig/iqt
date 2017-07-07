% COMPUTE_MODEL_DTI  A script to estimate DTI from DWIs.
%   This script is useful for computing DTIs on your own datasets. These
%   DTIs can then be super-resolved using RECONSTRUCT_HIRES_DTI.
%
%   This is a script, you have to edit the settings.
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
inp_dir = '/Data/root/'; % dir where DWI data is stored (eg your data root)
out_dir = '/Results/root/';  % typically root dir where results are stored
% list of test data subjects
data_folders = {'904044', '165840'}; %, '889579', '713239', '899885', '117324', '214423', '857263'};

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
dt_pref = 'dt_b1000_'; % DTI name prefix

% -- end of settings --

%%
open_matlabpool();


%% Estimation
compute_dti(inp_dir, out_dir, data_folders, sub_path, ...
            dw_file, bvals_file, bvecs_file, mask_file, ...
            dt_pref);
        
%%        
close_matlabpool();
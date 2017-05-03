% COMPUTE_MODEL_DTI  A script to estimate DTI from DWIs.
%   This script is useful for computing DTIs on your own datasets. These
%   DTIs can then be super-resolved using the IQT code.
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
addpath(genpath('.'));

% Set paths (always end directory paths with a forward/back slash)
inp_dir = ''; % dir where DWI data is stored (eg HCP data root)
out_dir = inp_dir; % dir where DTIs will be saved (default input dir)
% list of testing subjects
data_folders = {'904044', '165840'}; 

% Check
if strcmp(inp_dir, '')
    error('[IQT] Input DWI data root missing, please check paths in settings.');
end

% Optional settings
sub_path = 'T1w/Diffusion/'; % internal directory structure
dw_file = 'data.nii'; % DWI file
bvals_file = 'bvals'; % b-values file
bvecs_file = 'bvecs'; % b-vectors file
mask_file = 'nodif_brain_mask.nii'; % mask file
grad_file = ''; % gradient non-linearities (HCP only: grad_dev.nii)
                            % For non-HCP: grad_file = ''
dt_pref = 'dt_b1000_'; % DTI name prefix


%% Estimation
compute_dti(inp_dir, out_dir, data_folders, sub_path, ...
            dw_file, bvals_file, bvecs_file, mask_file, grad_file, ...
            dt_pref);
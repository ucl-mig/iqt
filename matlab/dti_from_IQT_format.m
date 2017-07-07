% DTI_FROM_IQT_FORMAT  A script to convert the super-resolved DTI
%   from IQT compatible format to 4D NIFTI images.
%   Default element ordering is MRtrix compatible.
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
inp_dir = '/Results/root/'; % dir where super-resolved DTI data is stored
out_dir = inp_dir;  % typically root dir where output will be stored
% list of test data subjects
data_folders = {'904044', '165840'}; %, '889579', '713239', '899885', '117324', '214423', '857263'};

% Check
if strcmp(inp_dir, '') || strcmp(out_dir, '')
    error('[IQT] Input/Output root missing, please check paths in settings.');
end

% Optional settings
sub_path = 'T1w/Diffusion/'; % internal directory structure
dt_pref = 'dt_b1000_recon_'; % input, IQT compatible format DTI name prefix
dti_file = 'dti_hires_b1000.nii'; % output DTI file (4D NIFTI)
% DTI element order
dti_ord = [1, 4, 6, 2, 3, 5]; %MRtrix compatible order (Dxx, Dyy, Dzz, Dxy, Dxz, Dyz)
%dti_ord = [1, 2, 3, 4, 5, 6]; %DTI element order (Dxx, Dxy, Dxz, Dyy, Dyz, Dzz)

% -- end of settings --


%% Convert
check_path(inp_dir);
check_path(out_dir);
check_path(sub_path);

for fi = 1:length(data_folders)
    if(~exist([out_dir data_folders{fi}  '/' sub_path], 'dir'))
        mkdir(out_dir, [data_folders{fi}  '/' sub_path]);
    end
    
    input_folder = [inp_dir data_folders{fi} '/' sub_path];
    output_folder = [out_dir data_folders{fi} '/' sub_path];
    
    % Read in the DTI data
    fprintf('Loading DTI: %s (%s)\n', dt_pref, data_folders{fi});
    [dt, hdr] = ReadDT_Volume([ input_folder dt_pref]);
    hdr.dime.dim(1) = 4;
    hdr.dime.dim(5) = 6;
    
    %reorder
    dt = dt(:,:,:,dti_ord+2);
    
    fprintf('Converting DTI subject: %s\n', data_folders{fi});
    nii = make_nii(dt);
    if ~isfield(hdr.hist, 'originator') && isfield(nii.hdr.hist, 'originator')
        hdr.hist.originator = nii.hdr.hist.originator;
    end
    nii.hdr = hdr;
    save_nii(nii, [output_folder dti_file]);

end

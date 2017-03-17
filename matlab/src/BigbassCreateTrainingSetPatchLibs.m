% Adapted from script_create_feature_vectors on amalga.
% Diverse training set.
data_folders = {'992774', '125525', '205119', '133928', '570243', '448347', '654754', '153025'};
%OutputPath = '/cs/research/vision/hcp/DCA_HCP.2013.3_Proc';
OutputPath = '/SAN/vision/hcp/DCA_HCP.2013.3_Proc';
%SubPath = 'T1w/Diffusion/';
SubPath = 'Diffusion/Diffusion/';
%SubPath = '';
dt_name = 'dt_b1000_' ;
%dt_name = 'dt_b2973_' ;
%dt_name = 'dt_all_' ;
%h4_name = 'h4_all_' ;

addpath('~/MSR_SuperRes');
addpath('~/MSR_SuperRes/nifti');

parfor fi = 1:length(data_folders)
    exec_BigbassCreateTrainingSetPatchLibs(fi, OutputPath, data_folders, SubPath, dt_name);
%   exec_BigbassCreateTrainingSetPatchLibsH4_2_H4(fi, OutputPath, data_folders, SubPath, dt_name, h4_name);
    
end



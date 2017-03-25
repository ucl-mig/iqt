% Testing on the cluster:

% %% Preprocessing:
% % Set the paths:
% dwi_dir = '/SAN/vision/hcp/HCP'; 
% dti_dir = '/SAN/vision/hcp/IQT-pre-release/DTI';
% traindata_dir = '/SAN/vision/hcp/IQT-pre-release/train-data';
% 
% settings.dt_name = 'dt_b1000_';
% settings.input_radius = 2; % the radius of the low-res input patch.
% settings.upsample_rate = 2; % the upsampling rate
% settings.subsample_rate = 32; % the rate of subsampling.
% settings.no_rnds = 8; % no of randomisations
% 
% 
% % Step 1: model computation (e.g. DTI) from DWIs.
% % compute the high-res and low-res DTI.
% data_folders = {'992774/T1w/Diffusion', '125525/T1w/Diffusion', '205119/T1w/Diffusion', '133928/T1w/Diffusion',...
%                 '570243/T1w/Diffusion', '448347/T1w/Diffusion', '654754/T1w/Diffusion', '153025/T1w/Diffusion',...
%                 '117324/T1w/Diffusion','904044/T1w/Diffusion'}; 
%             
% compute_dti(dwi_dir, dti_dir, data_folders, settings)
% 
% % Step 2: create a library of low-res and high-res patches.
% data_train_folders = {'992774/T1w/Diffusion', '125525/T1w/Diffusion', '205119/T1w/Diffusion', '133928/T1w/Diffusion',...
%                       '570243/T1w/Diffusion', '448347/T1w/Diffusion', '654754/T1w/Diffusion', '153025/T1w/Diffusion'};
% 
% compute_patchlib(dti_dir, dti_dir, data_train_folders, settings)
% 
% % Step 3: create a training set. 
% create_trainingset(dti_dir, traindata_dir, data_train_folders, settings)
% 

%% Train trees:

%Define parameters:

%Header of DT images.
settings.dt_name = 'dt_b1000_';
%Radius of the low-res input patch.
settings.input_radius = 2;
%Upsampling rate
settings.upsample_rate = 2; 
%Subsampling rate.
settings.subsample_rate = 32; 
%No of training sets. You train one tree on each set.
settings.no_rnds = 8; 
%Feature set used in Neuroimage paper. See PatchFeatureList.m for details.
settings.feature_version = 6; 


%Set the paths:

%Dir where training data is stored
traindata_dir = '/SAN/vision/hcp/IQT-pre-release/train-data'; 
%Dir where you save the trained trees
trees_dir = '/SAN/vision/hcp/IQT-pre-release/trees';

%Train trees
train_trees(traindata_dir, trees_dir, settings)

%% Testing (reconstruction):

%Header of DT images.
settings.dt_name = 'dt_b1000_';
%Radius of the low-res input patch.
settings.input_radius = 2;
%Upsampling rate
settings.upsample_rate = 2; 
%Subsampling rate.
settings.subsample_rate = 32; 
%Feature set used in Neuroimage paper. See PatchFeatureList.m for details.
settings.feature_version = 6; 
% Set true to perform boundary completion.
settings.edge = 0;

% Paths:
settings.input_dir = '/SAN/vision/hcp/IQT-pre-release/DTI'; 
settings.output_dir = '/SAN/vision/hcp/IQT-pre-release/recon';
settings.trees_dir = '/SAN/vision/hcp/IQT-pre-release/trees';
settings.trees_list = [1:8];
settings.patchlibs_dir = '/SAN/vision/hcp/IQT-pre-release/train-data';
data_folders = {'117324/T1w/Diffusion','904044/T1w/Diffusion'};

reconstruct_randomforests(data_folders, settings)
%visualise_results(data_folders, settings)



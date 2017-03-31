%% Train trees:

%Define parameters (default values below):
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
%Dir where training data (patch-libraries) are stored
traindata_dir = 'path/to/training/data'; 
%Dir where you save the trained trees
trees_dir = 'where/you/want/to/save/trees';

%Train trees
train_trees(traindata_dir, trees_dir, settings)
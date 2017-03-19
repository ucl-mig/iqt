% Train trees:
settings.dt_name = 'dt_b1000_';
settings.input_radius = 2; % the radius of the low-res input patch.
settings.upsample_rate = 2; % the upsampling rate
settings.subsample_rate = 32; % the rate of subsampling.
settings.no_rnds = 8; % no of training sets. You train one tree on each set.


% Set the paths:
traindata_dir = '~/tmp/iqt_codes/training_data';
trees_dir = '~/tmp/iqt_codes/trees';

train_trees(traindata_dir, trees_dir, settings)
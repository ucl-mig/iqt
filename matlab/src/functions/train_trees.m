function train_trees(train_dir, ds_rate, input_radius, ...
                     sample_rate, no_rnds, fv)
% TRAIN_TREES is actual function that trains random forest trees from the
%   data available in the train_dir directory.
%
%   Args:
%       TRAIN_DIR: This directory which contains all training outputs
%       DS_RATE: Super-resolution factor
%       INPUT_RADIUS: the input is a cubic patch of size (2*INPUT_RADIUS+1)^3
%       SAMPLE_RATE: data sampling rate for random sub-sampling
%       NO_RNDS: number of random sub-samples to create these many training
%                datasets
%       FV: feature version used for computing features for tree training
%
%   (always end directory paths with a forward/back slash)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

%% Define variables and display for checking
% Patch size (input patch radius). (!) Duplication. Remove
settings.n=input_radius; 
% Downsampling factor 
settings.ds=ds_rate;
% Output patch radius.
settings.m=settings.ds;
% Version of feature set
settings.fv=fv;
settings.sample_rate=sample_rate;
% Don't include spatial locations in the feature list.
settings.spatial=0;
% Type of tree truncation
settings.trunc = 'validation';
% Type of parameter mapping
settings.parmap = 'other';
%settings.parmap = 'H4_2_H4';
% The scale to set the value of the constant term in each linear
% transformation.
settings.scale = 1E-3;
if(strcmp(settings.parmap, 'H4_2_H4'))
    settings.scale = 1;
    settings.noRankTest = 1;
end
% Decision on whether to insist on full rank S in SplitInfoGain. Used this
% for all DTI patch training, but need to drop if for H4_2_H4. Probably
% somewhere in between is generally applicable.
settings.noRankTest = 0;
if(strcmp(settings.parmap, 'H4_2_H4'))
    settings.noRankTest = 1;
end
% Whether to suppress the warnings that occur due to inversion of low rank
% matrices, which rarely affect the outcome, as they usually occur for
% thresholds that are not eventually selected.
settings.suppress_low_rank_warnings = 1;
% The maximum number of data points to use during node threshold selection.
% For DTI.
settings.max_size_for_split = 100000;
% For H4.
if(strcmp(settings.parmap, 'H4_2_H4'))
    settings.max_size_for_split = 1000000;
end

disp(settings)


%% Get the list of training sets
TrainingDataFiles = {};
for si = 1:no_rnds
    TrainingDataFiles{end+1} = sprintf('_DS%02i_%ix%ix%i_%ix%ix%i_Sub%03i_%04i.mat', ...
    settings.ds, 2*settings.n+1,2*settings.n+1,2*settings.n+1, settings.m, settings.m, settings.m, settings.sample_rate, si);
end


%% Train trees.
% loop over training data sets.
for tfi=1:length(TrainingDataFiles)
    fprintf('Training tree %i/%i\n', tfi, length(TrainingDataFiles));
    fprintf('on dataset: PatchLibs%s\n', TrainingDataFiles{tfi});
    load([train_dir '/PatchLibs' TrainingDataFiles{tfi}])

    % Compute candidate split features.
    FeatureFileName = [train_dir '/FeaturesV' int2str(settings.fv) TrainingDataFiles{tfi}];
    if(~exist(FeatureFileName, 'file'))
        display('Computing features...')
        tic;
        if(strcmp(settings.parmap, 'H4_2_H4'))
            feature_source = comipatchlib(:,1:((2*settings.n+1)^3*6));
        else
            feature_source = comipatchlib;
        end
        features = PatchFeatureList(feature_source, settings.n, settings.ds, settings.fv, settings.spatial, patchlibindices);
        save(FeatureFileName, 'features');
        feature_source = [];
        toc
    else
        display('Loading features...')
        load(FeatureFileName);
    end

    % Train a regression tree with validation set
    display('Start training ...')
    train_inds = find(rand(1,length(features))<0.5);
    test_inds = setdiff(1:length(features), train_inds);

    train_features = features(train_inds,:);
    %train_in = comipatchlib(train_inds,:);
    if(strcmp(settings.parmap, 'H4_2_H4'))
        train_in = comipatchlib(train_inds,((2*settings.n+1)^3*6 + 1):end);
        settings.input_patch_scaling = mean(abs(train_in));
        train_in = train_in./repmat(settings.input_patch_scaling, [size(train_in,1), 1]);
        train_in = [train_in'; settings.scale*ones(1,length(train_inds))]';
        
        train_out = comopatchlib(train_inds,:);
        settings.output_patch_scaling = mean(abs(train_out));
        train_out = train_out./repmat(settings.output_patch_scaling, [size(train_out,1), 1]);
        
    else
        train_in = [comipatchlib(train_inds,:)'; settings.scale*ones(1,length(train_inds))]';
        train_out = comopatchlib(train_inds,:);
    end

    test_features = features(test_inds,:);
    %test_in = comipatchlib(test_inds,:);
    if(strcmp(settings.parmap, 'H4_2_H4'))
        test_in = comipatchlib(test_inds,((2*settings.n+1)^3*6 + 1):end);
        test_in = test_in./repmat(settings.input_patch_scaling, [size(test_in,1), 1]);
        test_in = [test_in'; settings.scale*ones(1,length(test_inds))]';
        
        test_out = comopatchlib(test_inds,:);
        test_out = test_out./repmat(settings.output_patch_scaling, [size(test_out,1), 1]);
    else
        test_in = [comipatchlib(test_inds,:)'; settings.scale*ones(1,length(test_inds))]';
        test_out = comopatchlib(test_inds,:);
    end
    
    % Compute the root node of the tree: perform linear regression
    M = train_in\train_out;
    S = train_out - (M'*train_in')';
    globalE = S;
    S = S'*S/length(train_in);
    Sv = test_out - (M'*test_in')';
    Sv = Sv'*Sv/length(test_in);
    
    tic;
    rootnode = struct('IsLeaf', 1,...
                      'SplitFeatureIndex', 0,...
                      'FeatureThreshold', 0,...
                      'LeftChildIndex', 0,...
                      'RightChildIndex', 0,...
                      'Transformation', struct('M', M, 'S', S, 'Sv', Sv, 'RegressionType', 'linear'),...
                      'InfoContent', 2.0*sum(log(diag(chol(S))))*length(train_features),...
                      'TrainingDataIndices', 1:length(train_features),...
                      'Visited', 0);

    if(strcmp(settings.trunc, 'validation'))
        rootnode.TestDataIndices=1:length(test_features);
    end
    
    % Grow the tree:
    tree = {rootnode};
    nextind = 2;
    max_layers = 20;
    for i=1:max_layers
        display(sprintf('Layer %i', i));
        change = 0;
        for j=1:length(tree)
            if(tree{j}.IsLeaf && ~tree{j}.Visited)
                display(sprintf('Processing node %i (%i data points)', j, length(tree{j}.TrainingDataIndices)));

                if(strcmp(settings.trunc, 'validation'))
                    % Validation set splitting
                    [tn, ln, rn, split_rel_error] = SplitNodeValidPar(tree{j}, train_in, train_out, train_features, test_in, test_out, test_features, settings);
                end
                
                display(sprintf('Relative split error: %d.', split_rel_error));
                
                info_gain = tn.InfoContent-ln.InfoContent-rn.InfoContent;
                if(info_gain>0 && split_rel_error<0)

                    tree{j} = tn;
                    tree{j}.IsLeaf = 0;
                    tree{j}.LeftChildIndex = nextind;
                    tree{j}.RightChildIndex = nextind+1;
                    tree{nextind} = ln;
                    tree{nextind+1} = rn;
                    nextind = nextind+2;
                    change = 1;
                end
                tree{j}.Visited = 1;
            end
        end
        
        if(~change)
            break;
        end
    end
    toc
    
    % Save the tree:
    if ~exist(train_dir, 'dir')
        mkdir(train_dir);
    end
    display('Saving the tree ...')
    tree = strip_tree(tree);
    save([train_dir '/RegTreeValV' int2str(settings.fv) TrainingDataFiles{tfi}], 'tree', 'settings');
    %save([output_dir '/RegTreeValV' int2str(settings.fv) TrainingDataFiles{tfi}], 'tree', 'train_inds', 'settings');
    fprintf('\n')
end

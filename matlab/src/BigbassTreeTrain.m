
%% Load in training data
DATAPATH='/cs/research/vision/hcp/DCA_HCP.2013.3_Proc/TrainingData';
DATAPATH='/SAN/vision/hcp/DCA_HCP.2013.3_Proc/TrainingData';


%% Training sets for generalisation demonstration
% Super-res runs                
TrainingDataFiles = {'DiverseDS02_5x5_2x2_TS8_SRi032_0001.mat', 'DiverseDS02_5x5_2x2_TS8_SRi032_0002.mat', 'DiverseDS02_5x5_2x2_TS8_SRi032_0003.mat', 'DiverseDS02_5x5_2x2_TS8_SRi032_0004.mat',...
                     'DiverseDS02_5x5_2x2_TS8_SRi032_0005.mat', 'DiverseDS02_5x5_2x2_TS8_SRi032_0006.mat', 'DiverseDS02_5x5_2x2_TS8_SRi032_0007.mat', 'DiverseDS02_5x5_2x2_TS8_SRi032_0008.mat'};

%% Define variables and display for checking
% Patch size
settings.n=3;
% Downsampling factor
settings.ds=2;
% Feature version
settings.fv=6;
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

settings
    

%% Loop over training data sets.
for tfi=1:length(TrainingDataFiles)

    display(TrainingDataFiles{tfi})
    load([DATAPATH '/PatchLibs' TrainingDataFiles{tfi}])

    % Compute candidate split features.
    FeatureFileName = [DATAPATH '/FeaturesV' int2str(settings.fv) '_' TrainingDataFiles{tfi} '.mat'];
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
        load(FeatureFileName);
    end

    % Train a regression tree with validation set
    display('Training tree...')
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

    M = train_in\train_out;
    S = train_out - (M'*train_in')';
    % Maintain a global error score for AIC/BIC implementation
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
                else
                    % AIC or BIC splitting
                    [tn, ln, rn, newGlobalE] = SplitNodeAIC_Par(tree{j}, train_in, train_out, train_features, globalE);
                    SigNoSplit = sum(sum(globalE.^2))/length(train_inds)
                    SigAfterSplit = sum(sum(newGlobalE.^2))/length(train_inds)
                    if(strcmp(settings.trunc, 'AIC'))
                        delAIC = length(train_inds)*(log(SigAfterSplit) - log(SigNoSplit)) + 2*(settings.ds*settings.ds*settings.ds+2)*(2*settings.n+1)^3 + 4;
                        split_rel_error = delAIC;
                    elseif(strcmp(settings.trunc, 'BIC'))
                        delBIC = length(train_inds)*(log(SigAfterSplit) - log(SigNoSplit)) + log(length(train_inds))*((settings.ds*settings.ds*settings.ds+2)*(2*settings.n+1)^3 + 2);
                        split_rel_error = delBIC;
                    else
                        error('Unknown truncation strategy');
                    end
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
                    if(strcmp(settings.trunc, 'AIC') || strcmp(settings.trunc, 'BIC'))
                        globalE = newGlobalE;
                    end
                    
                end
                tree{j}.Visited = 1;
            end
        end
        %save /tmp/tree_temp.mat tree
        if(~change)
            break;
        end
    end

    toc

    save([DATAPATH '/RegTreeValV' int2str(settings.fv) '_' TrainingDataFiles{tfi}], 'tree', 'train_inds', 'settings');

end

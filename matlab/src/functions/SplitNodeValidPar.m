function [pn, ln, rn, split_rel_error] = SplitNodeValidPar(parent_node, train_in, train_out, train_features, test_in, test_out, test_features, settings, features_list)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

if (nargin < 9)
    features_list = 1:size(train_features,2); %optimise over all features unless specified! 
end

pn = parent_node;
parfor i = 1:length(features_list) %Optimize the splitting for each feature. 
    display(sprintf('     Computing optimal splitting over %d th feature',features_list(i))) 
    [pnt, lnt, rnt] = SplitNodeOptGivenFeature(parent_node, train_in, train_out, train_features, features_list(i), settings);
    info_gain = pnt.InfoContent-lnt.InfoContent-rnt.InfoContent; 
    display(sprintf('     with information gain %d',info_gain)) 
    
    pns{i} = pnt;
    lns{i} = lnt;
    rns{i} = rnt;
    igs(i) = info_gain;
end
%igs;
[max_info_gain, max_ig_ind] = max(igs); %select the feature which maximizes the information gain.
pn = pns{max_ig_ind};
ln = lns{max_ig_ind};
rn = rns{max_ig_ind};


% Test performance on the validation dataset.
split_rel_error = 1;

if(max_info_gain > 0)
    test_out_pn = NodeRegress(pn, test_in(pn.TestDataIndices,:));
    error_parent = sum(sum((test_out(pn.TestDataIndices,:) - test_out_pn).^2));

    test_features_pn = test_features(pn.TestDataIndices,:);

    test_subset_left_inds = find(test_features_pn(:,pn.SplitFeatureIndex)<pn.FeatureThreshold);
    test_subset_right_inds = setdiff(1:size(test_features_pn,1),test_subset_left_inds);

    ln.TestDataIndices = pn.TestDataIndices(test_subset_left_inds);
    rn.TestDataIndices = pn.TestDataIndices(test_subset_right_inds);

    test_out_left = NodeRegress(ln, test_in(ln.TestDataIndices,:));
    test_out_right = NodeRegress(rn, test_in(rn.TestDataIndices,:));
    
    % Add in the error covariances from the validation set
    Sv_left = test_out(ln.TestDataIndices,:) - test_out_left;
    ln.Transformation.Sv = Sv_left'*Sv_left/length(ln.TestDataIndices);
    Sv_right = test_out(rn.TestDataIndices,:) - test_out_right;
    rn.Transformation.Sv = Sv_right'*Sv_right/length(rn.TestDataIndices);

    test_out_split = zeros(size(test_out(pn.TestDataIndices,:)));
    test_out_split(test_subset_left_inds,:) = test_out_left;
    test_out_split(test_subset_right_inds,:) = test_out_right;

    error_split = sum(sum((test_out(pn.TestDataIndices,:) - test_out_split).^2));

    split_rel_error = error_split - error_parent;
end



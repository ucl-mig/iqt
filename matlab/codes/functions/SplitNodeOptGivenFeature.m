function [pn, ln, rn] = SplitNodeOptGivenFeature(parent_node, train_in, train_out, train_features, feature_ind, settings)

% Don't use all the training data for finding the split if the training set
% it large.  pn is a temporary version of parent_node with a subsampled set
% of training data.
pn = parent_node;
if(length(pn.TrainingDataIndices)>settings.max_size_for_split)
    [a b] = sort(rand(1,length(pn.TrainingDataIndices)));
    pn.TrainingDataIndices = pn.TrainingDataIndices(b(1:settings.max_size_for_split));
    display(sprintf('Using subset of %i of %i data points for feature selection.', settings.max_size_for_split, length(pn.TrainingDataIndices)));
end

% Bounds on range are the smallest and largest values appearing in the
% list.
upperbound = max(train_features(pn.TrainingDataIndices,feature_ind));
lowerbound = min(train_features(pn.TrainingDataIndices,feature_ind));

% Now search for the split threshold using the subsampled training data
% set.
[max_thresh, neg_info_gain, exitflag, output] = fminbnd(@(thresh) SplitInfoGain(thresh, train_features, feature_ind, pn, train_in, train_out, settings), lowerbound, upperbound);

% Now compute the information gain and transformation using the full
% training data available to this node.
[neg_info_gain, M_left, M_right, S_left, S_right, train_inds_left, train_inds_right] = SplitInfoGain(max_thresh, train_features, feature_ind, parent_node, train_in, train_out, settings);
info_gain = -neg_info_gain;

pn = parent_node;
ln = struct('InfoContent', parent_node.InfoContent/2);
rn = struct('InfoContent', parent_node.InfoContent/2);

if(info_gain>0)    
    ln = parent_node;
    ln.Transformation.M = M_left;
    ln.Transformation.S = S_left;
%   The first estimate of the log det breaks down for large matrices.
%   ln.InfoContent = log(det(S_left*iscale))*length(train_inds_left);
    ln.InfoContent = length(train_inds_left)*2.0*sum(log(diag(chol(S_left))));
    ln.TrainingDataIndices = train_inds_left;
    
    rn = parent_node;
    rn.Transformation.M = M_right;
    rn.Transformation.S = S_right;
%   rn.InfoContent = log(det(S_right*iscale))*length(train_inds_right);
    rn.InfoContent = length(train_inds_right)*2.0*sum(log(diag(chol(S_right))));
    rn.TrainingDataIndices = train_inds_right;
    
    pn.IsLeaf = 0;
    pn.SplitFeatureIndex = feature_ind;
    pn.FeatureThreshold = max_thresh;
end
    


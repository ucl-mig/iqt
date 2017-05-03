% Computes the information gain from splitting a node on a particular
% feature at a particular threshold and constructing linear mappings for
% each subset instead of a linear mapping for all the data.
function [neg_info_gain, M_left, M_right, S_left, S_right, train_inds_left, train_inds_right] = SplitInfoGain(thresh, train_features, feature_ind, parent_node, train_in, train_out, settings)


% Minimum size of data set for fitting a linear model.
min_datasize = size(parent_node.Transformation.M,1);
% For constant, we need at least 2 points to compute the information
% content.


% Divide into subsets based on threshold
subset_left_inds = find(train_features(parent_node.TrainingDataIndices,feature_ind)<thresh);
subset_right_inds = setdiff(1:length(parent_node.TrainingDataIndices),subset_left_inds);

% Default assignments
% Initialize the information gain to something smaller than it could ever
% be so that the search doesn't get stuck in unviable regions of the
% threshold parameter.
info_gain = -sum(sum(train_in))^2*1E4;
M_left = 0;
M_right = 0;
S_left = 0;
S_right = 0;
train_inds_left = 0;
train_inds_right = 0;

% Compute transformations and information gain if possible
if(length(subset_left_inds)>min_datasize && length(subset_right_inds)>min_datasize)
    train_inds_left = parent_node.TrainingDataIndices(subset_left_inds);
    train_in_left  = train_in(train_inds_left, :);
    train_out_left = train_out(train_inds_left, :);
  
    if(settings.suppress_low_rank_warnings)
            warning('off', 'MATLAB:rankDeficientMatrix');
    end
    
    M_left = train_in_left\train_out_left;%pinv(train_in_left)*train_out_left;
    S_left = train_out_left - (M_left'*train_in_left')';

%train_error = sum(sum(S_left.^2));
    S_left = S_left'*S_left/length(train_inds_left);

    train_inds_right = parent_node.TrainingDataIndices(subset_right_inds);
    train_in_right = train_in(train_inds_right, :);
    train_out_right = train_out(train_inds_right, :);
    
    if(settings.suppress_low_rank_warnings)
        warning('off', 'MATLAB:rankDeficientMatrix');
    end
    
    M_right = train_in_right\train_out_right;%pinv(train_in_right)*train_out_right;
    S_right = train_out_right - (M_right'*train_in_right')';
    
%train_error = train_error + sum(sum(S_right.^2));
    S_right = S_right'*S_right/length(train_inds_right);

    % Note that if the number of training points is reduced through
    % randomization, as it is during threshold identification for large
    % training sets, this is not strictly the information gain, because the
    % value stored in the parent_node assumes a different number of total
    % training points.  However, the difference is a constant so the outcome
    % of the threshold selection is unaffected.
    %info_parent = length(train_in)*2*sum(log(diag(chol(parent_node.Transformation.S))))
    if(settings.noRankTest || (rank(S_left)==length(S_left) && rank(S_right)==length(S_right)))
        try 
            info_left = length(train_in_left)*2.0*sum(log(diag(chol(S_left))));
            info_right = length(train_in_right)*2.0*sum(log(diag(chol(S_right))));
            info_gain = parent_node.InfoContent-info_left-info_right;
        catch err
            display([err.message ' In SplitInfoGain use of chol. Feature ' int2str(feature_ind) '; threshold ' num2str(thresh) '.'])
            info_gain = -1;
        end

%[thresh length(train_in_left) length(train_in_right) info_left info_right info_gain train_error]
    end

end

neg_info_gain = -info_gain;


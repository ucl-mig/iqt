% Performs regression for an input set of data points using a single tree
% node.  The form of regression depends on the contents of the node so may
% be, for example, linear or constant.
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

function out = NodeRegress(node, in)

if(~isfield(node.Transformation, 'RegressionType'))
    node.Transformation.RegressionType = 'linear';
end

if(strcmp(node.Transformation.RegressionType, 'linear'))
    out = (node.Transformation.M'*in')';
elseif(strcmp(node.Transformation.RegressionType, 'constant'))
    out = repmat(node.Transformation.M, [size(in, 1), 1]);
else
    error(['Unknown regression type: ' node.Transformation.RegressionType '.']);
end



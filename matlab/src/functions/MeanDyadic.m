function Y = MeanDyadic(pds)
% Computes the dyadic tensor of a list of directions to enable
% a measure of direction concentration.
%
% The dyadic tensor is Y = 1/N \sum_i^N n_i n_i^T where n_i is the
% i-th of N directions.
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

Y = zeros(3,3);
for i=1:length(pds)
  Y = Y + pds(i,:)'*pds(i,:);
end

Y = Y/length(pds);


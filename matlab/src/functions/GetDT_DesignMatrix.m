function X = GetDT_DesignMatrix(bvals, bvecs, L)

% If no gradient non-linearity correction provided, just use the unaltered
% b matrix.
if(nargin<3 || numel(L)<=1)
    new_bvecs = bvecs;
    new_bvals = bvals;
else
    I = eye(3);

    % correct bvecs and calculate their norm
    v = (I+L)*bvecs;
    nrm = sqrt(sum(v.^2));

    % rotate bvecs to voxel coord system 
    % and correct both bvecs and bvals
    new_bvecs = v./repmat(nrm,3,1);
    new_bvals = nrm.^2.*bvals;
end

% Construct the design matrix for linear estimation of the DT.
X = [ones(1, length(new_bvals));...
    -new_bvals.*new_bvecs(1,:).*new_bvecs(1,:);...
    -2*new_bvals.*new_bvecs(1,:).*new_bvecs(2,:);...
    -2*new_bvals.*new_bvecs(1,:).*new_bvecs(3,:);...
    -new_bvals.*new_bvecs(2,:).*new_bvecs(2,:);...
    -2*new_bvals.*new_bvecs(2,:).*new_bvecs(3,:);...
    -new_bvals.*new_bvecs(3,:).*new_bvecs(3,:)]';



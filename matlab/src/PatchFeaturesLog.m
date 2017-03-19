% Computes various features for each tensor in a tensor patch.  This
% version operates on the log tensors, but computes the same set of
% features as PatchFeatures - only the input is different.  The
% features are:
%
% - Direction concentration of the tensor principal directions.
% - Average principle direction.
% - Mean first eigenvalue.
% - Std first eigenvalue.
% - Magnitude of first eigenvalue (L1) gradient over patch.
% - Absolute dot product of L1 gradient orientation with average 
% principle direction.
% - Then the four features above for E1 repeated for each of the
% following:
%   - Second eigenvalue E2.
%   - Third eigenvalue E3.
%   - Tensor linearity (E1-E2)/E1.
%   - Tensor planarity (E2-E3)/E1.
%   - Tensor isotropy E3/E1.
function fts = PatchFeaturesLog(patch)

p=3;
rpatch = reshape(patch, [p,p,p,6]);
for i=1:p
    for j=1:p
        for k=1:p
            ldt = MakeDT_Matrix(rpatch(i,j,k,1),rpatch(i,j,k,2),rpatch(i,j,k,3),rpatch(i,j,k,4),rpatch(i,j,k,5),rpatch(i,j,k,6));
            [R E] = eig(ldt);
            E = diag(exp(diag(E)));
            % Eigensystem
            esys(i,j,k,1:6) = [diag(E)' R(:,3)'];
            % Linearity
            esys(i,j,k,7) = (E(3,3)-E(2,2))/E(3,3);
            % Planarity
            esys(i,j,k,8) = (E(2,2)-E(1,1))/E(3,3);
            % Isotropy
            esys(i,j,k,9) = E(1,1)/E(3,3);
        end
    end
end

% Include the mean of each scalar feature
meanf = ones(3,3,3)/27;

% 3D gradient operators
sobel3x(:,:,1) = [[1 0 -1];[2 0 -2];[1 0 -1]];
sobel3x(:,:,2) = [[2 0 -2];[4 0 -4];[2 0 -2]];...
sobel3x(:,:,3) = [[1 0 -1];[2 0 -2];[1 0 -1]];

sobel3y(:,:,1) = sobel3x(:,:,1)';
sobel3y(:,:,2) = sobel3x(:,:,2)';
sobel3y(:,:,3) = sobel3x(:,:,3)';

sobel3z(:,:,1) = sobel3x(:,1,:);
sobel3z(:,:,2) = sobel3x(:,2,:);
sobel3z(:,:,3) = sobel3x(:,3,:);

% First compute statistics of the set of principle directions.
ind = 1;
pds = reshape(esys(:,:,:,4:6),[27,3]);
Y = MeanDyadic(pds);
[R E] = eig(Y);
meanPD = R(:,3);
dirconc = E(3,3);
fts(ind) = dirconc;
ind = ind + 1;
fts(ind:(ind+2)) = meanPD;
ind = ind + 3;

% Now statistics of the various orientationally invariant indices.
for i=[1 2 3 7 8 9]
    
    % Mean of index.
    fts(ind) = sum(sum(sum(meanf.*esys(:,:,:,i))));
    ind = ind + 1;
    
    % Standard deviation of index.
    fp = squeeze(esys(:,:,:,i));
    fts(ind) = std(fp(:));
    ind = ind + 1;

    % Index gradient over patch.
    sx = sum(sum(sum(sobel3x.*esys(:,:,:,i))));
    sy = sum(sum(sum(sobel3y.*esys(:,:,:,i))));
    sz = sum(sum(sum(sobel3z.*esys(:,:,:,i))));

    gmag = sqrt(sx^2 + sy^2 + sz^2);
    fts(ind) = gmag;
    ind = ind + 1;
    
    % Gradient direction relative to mean principal direction
    graddir = [sx, sy, sz]/gmag;
    agdotd = abs(sum(graddir'.*meanPD));
    fts(ind) = agdotd;
    ind = ind + 1;
end


    
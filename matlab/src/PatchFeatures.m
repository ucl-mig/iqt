% Computes various features for each tensor in a tensor patch.  The
% features are:
%
% 1 Direction concentration of the tensor principal directions.
% 2-4 Average principle direction (x,y,z)
% 5 Mean first eigenvalue (L1).
% 6 Std L1.
% 7 Magnitude of L1 gradient over patch.
% 8 Absolute dot product of L1 gradient orientation with average 
% principle direction.
% Then the four features above for L1 repeated for each of the
% following:
%   9-12 Second eigenvalue L2.
%   13-16 Third eigenvalue L3.
%   17-20 Tensor linearity (L1-L2)/L1.
%   21-24 Tensor planarity (L2-L3)/L1.
%   25-28 Tensor isotropy L3/L1.
function fts = PatchFeatures(patch)

ps=floor((length(patch)/6)^(1/3)+0.001);
rpatch = reshape(patch, [ps,ps,ps,6]);
% For these particular features, we just use the inner pxpxp block.
p = 3;
rpatch = rpatch(((ps-p)/2+1):((ps+p)/2),((ps-p)/2+1):((ps+p)/2),((ps-p)/2+1):((ps+p)/2),:);

for i=1:p
    for j=1:p
        for k=1:p
            ldt = MakeDT_Matrix(rpatch(i,j,k,1),rpatch(i,j,k,2),rpatch(i,j,k,3),rpatch(i,j,k,4),rpatch(i,j,k,5),rpatch(i,j,k,6));
            [R E] = eig(ldt);
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


    
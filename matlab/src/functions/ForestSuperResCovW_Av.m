function dtreconrf = ForestSuperResCovW_Av(dt_lowres, trees, n, m, patch_feature_map, us, rescale_factor)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

if(nargin<7)
    rescale_factor = 1;
end

for i=1:length(trees)
    
    % Do the reconstruction
    [dtreconrt, accum] = RegressionTreeSuperResOverlap(dt_lowres, trees{i}, n, m, patch_feature_map, us, 'cov_weighted_average', rescale_factor);

%     % Reweight the reconstruction by the accumulator for this contribution.
%     fg_inds = find(dtreconrt(:,:,:,1)==0);
%     for j=2:8;
%         v = dtreconrt(:,:,:,j);
%         v(fg_inds) = v(fg_inds).*accum(fg_inds);
%         dtreconrt(:,:,:,j) = v;
%     end
% 
    if(i==1)
        dtreconrf = dtreconrt;
        forest_accum = accum;
    else
        dtreconrf = dtreconrf+dtreconrt;
        forest_accum = forest_accum + accum;
    end

end

% Normalise by the summed covariances.
[RXSIZE, RYSIZE, RZSIZE, junk] = size(dtreconrf);
for i=1:RXSIZE
    for j=1:RYSIZE
        for k=1:RZSIZE
            if(dtreconrf(i,j,k,1)==0);
                dtreconrf(i,j,k,3:8) = squeeze(forest_accum(i,j,k,:,:))\squeeze(dtreconrf(i,j,k,3:8));
            end
        end
    end
end    


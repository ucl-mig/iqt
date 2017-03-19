function dtreconrf = ForestSuperResEdgeCovW_Av(dt_lowres, trees, n, m, patch_feature_map, us, T, comipatchmean)

for i=1:length(trees)
    
    % Do the reconstruction
    [dtreconrt, accum] = RegressionTreeSuperResOverlapEdge(dt_lowres, trees{i}, n, m, patch_feature_map, us, T, comipatchmean, 'cov_weighted_average');

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


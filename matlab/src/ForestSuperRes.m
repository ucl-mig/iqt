function dtreconrf = ForestSuperRes(dt_lowres, trees, n, m, patch_feature_map, us, rescale_factor)

if(nargin<7)
    rescale_factor = 1;
end

for i=1:length(trees)
    
    % Do the reconstruction
    [dtreconrt, accum] = RegressionTreeSuperResOverlap(dt_lowres, trees{i}, n, m, patch_feature_map, us, 'weighted_average', rescale_factor);

    % Reweight the reconstruction by the accumulator for this contribution.
    fg_inds = find(dtreconrt(:,:,:,1)==0);
    for j=2:8;
        v = dtreconrt(:,:,:,j);
        v(fg_inds) = v(fg_inds).*accum(fg_inds);
        dtreconrt(:,:,:,j) = v;
    end

    if(i==1)
        dtreconrf = dtreconrt;
        forest_accum = accum;
    else
        dtreconrf = dtreconrf+dtreconrt;
        forest_accum = forest_accum + accum;
    end

end

fg_inds = find(dtreconrf(:,:,:,1)==0);
for j=2:8;
    v = dtreconrf(:,:,:,j);
    v(fg_inds) = v(fg_inds)./forest_accum(fg_inds);
    dtreconrf(:,:,:,j) = v;
end


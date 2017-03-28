function [dtreconrf, forest_accum] = ForestSuperResEdge(dt_lowres, trees, n, m, patch_feature_map, us, T, comipatchmean, rescale_factor, scale_const, overlap)
% dtreconrf: estimated high-res dti
% forest_accum: average precision of the prediction over all trees
for i=1:length(trees)
     
    % Do the reconstruction
    [dtreconrt, accum] = TreeSuperResEdge(dt_lowres, trees{i}, n, m, patch_feature_map, us, T, comipatchmean, overlap, rescale_factor, scale_const);

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

%forest_accum = forest_accum/length(trees);

fg_inds = find(dtreconrf(:,:,:,1)==0);
for j=2:8;
    v = dtreconrf(:,:,:,j);
    v(fg_inds) = v(fg_inds)./forest_accum(fg_inds);
    dtreconrf(:,:,:,j) = v;
end


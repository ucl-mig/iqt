function [dtreconrf, mx_hash] = ForestSuperResEdge(dt_lowres, trees, n, m, patch_feature_map, us, T, comipatchmean, rescale_factor, scale_const, tree_settings, mx_hash)

if(nargin<12)
    mx_hash.n = n;
    mx_hash.m = m;
    mx_hash.T = T;
    mx_hash.comipatchmean = comipatchmean;
end

for i=1:length(trees)
    
    if(~isfield(tree_settings{i}, 'input_patch_scaling'))
        tree_settings{i}.input_patch_scaling = 1;
    end
    if(~isfield(tree_settings{i}, 'output_patch_scaling'))
        tree_settings{i}.output_patch_scaling = 1;
    end
    
    % Do the reconstruction
    [dtreconrt, accum, leafmap, mx_hash] = RegressionTreeSuperResOverlapEdge(dt_lowres, trees{i}, n, m, patch_feature_map, us, T, comipatchmean, 'weighted_average', rescale_factor, scale_const, tree_settings{i}.input_patch_scaling, tree_settings{i}.output_patch_scaling, mx_hash);

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


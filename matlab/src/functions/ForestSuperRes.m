function [dtreconrf, forest_accum, forest_var, leafmap] = ForestSuperRes(dt_lowres, trees, n, m, patch_feature_map, us, rescale_factor, ValidVar,  max_layer)

if(nargin<7)
    rescale_factor = 1;
end

if (nargin<9)
    max_layer = inf;
end

for i=1:length(trees)    
    % Do the reconstruction
    if ValidVar == 1
    [dtreconrt, accum, accum_var, leafmap] = TreeSuperRes(dt_lowres, trees{i}, n, m, patch_feature_map, us, 'weighted_average_valid', rescale_factor, max_layer);
    else
    [dtreconrt, accum, accum_var, leafmap] = TreeSuperRes(dt_lowres, trees{i}, n, m, patch_feature_map, us, 'weighted_average', rescale_factor, max_layer);
    end
    
    
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
        forest_var   = accum_var;
    else
        dtreconrf = dtreconrf+dtreconrt;
        forest_accum = forest_accum + accum;
        forest_var   = forest_var + accum_var;
    end
end

fg_inds = find(dtreconrf(:,:,:,1)==0); %average over trees!
    for j=2:8;
        v = dtreconrf(:,:,:,j);
        v(fg_inds) = v(fg_inds)./forest_accum(fg_inds);
        dtreconrf(:,:,:,j) = v;
    end



function fts = PatchFeaturesV6(patch)
% These are the PatchFeaturesV5 without the orientation-specific features.
% 1:3 Eigenvalues of DT at patch centre
% 4:7 Linearity, planarity, isotropy, trace.
% If n>0
% 8:10 Mean eigenvalues over central 3x3x3 patch
% 11 Orientation dispersion over central 3x3x3 patch
% 12:15 Mean linearity, planarity, isotropy, trace over central 3x3x3 patch.
% If n>1
% 16:18 Mean eigenvalues over full patch
% 19 Orientation dispersion over full patch
% 20:23 Mean linearity, planarity, isotropy, trace over full patch.
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

fts5 = PatchFeaturesV5(patch);
if(length(fts5)==10)
    fts = fts5([1:3,7:10]);
elseif(length(fts5)==21)
    fts = fts5([1:3,7:13,17:21]);
else
    fts = fts5([1:3,7:13,17:24,28:32]);
end



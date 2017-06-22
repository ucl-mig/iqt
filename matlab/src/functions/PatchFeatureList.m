function features = PatchFeatureList(patchlib, n, ds, fv, spatial, patchlibindices)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

if(fv==1)
    features = zeros(length(patchlib), 28);
    parfor i=1:length(patchlib)
        features(i,:) = PatchFeatures(patchlib(i,:));
    end
    
    % If the input patch size is 3x3, the last 18 features are the same as
    % features 7:24, so remove them to save computation.
    if(n==1)
        features = features(:,1:24);
    end
    
elseif(fv==2)
    features = zeros(length(patchlib), 42);
    parfor i=1:length(patchlib)
        features(i,:) = PatchFeaturesV2(patchlib(i,:));
    end
    
    % If the input patch size is 3x3, the last 18 features are the same as
    % features 7:24, so remove them to save computation.
    if(n==1)
        features = features(:,1:24);
    end

elseif(fv==4)
    if(n==0);
        features = zeros(length(patchlib), 6);
    elseif(n==1);
        features = zeros(length(patchlib), 13);
    else
        features = zeros(length(patchlib), 20);
    end
    parfor i=1:length(patchlib)
        features(i,:) = PatchFeaturesV4(patchlib(i,:));
    end
    
elseif(fv==5)
    if(n==0);
        features = zeros(length(patchlib), 10);
    elseif(n==1);
        features = zeros(length(patchlib), 21);
    else
        features = zeros(length(patchlib), 32);
    end
    parfor i=1:length(patchlib)
        features(i,:) = PatchFeaturesV5(patchlib(i,:));
    end
    
elseif(fv==6)
    if(n==0);
        features = zeros(length(patchlib), 7);
    elseif(n==1);
        features = zeros(length(patchlib), 15);
    else
        features = zeros(length(patchlib), 23);
    end
    parfor i=1:length(patchlib)
        features(i,:) = PatchFeaturesV6(patchlib(i,:));
    end
end

% Add spatial locations if specified
if(spatial)
    locations = PatchLibSpatialLocations(patchlibindices, ds, n);
    features = [features'; locations']';
end


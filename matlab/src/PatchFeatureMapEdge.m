function pfmap = PatchFeatureMapEdge(dt_lowres, n, T, comipatchmean, version)

if(nargin<4)
    version = 1;
end

if(version==1)
    numFeatures = 28;
elseif(version==2)
    if(n==1)
        numFeatures = 24;
    else
        numFeatures = 42;
    end
elseif(version==4)
    if(n==0)
        numFeatures = 6;
    elseif(n==1)
        numFeatures = 13;
    else
        numFeatures = 20;
    end
elseif(version==5)
    if(n==0)
        numFeatures = 10;
    elseif(n==1)
        numFeatures = 21;
    else
        numFeatures = 32;
    end
elseif(version==6)
    if(n==0)
        numFeatures = 7;
    elseif(n==1)
        numFeatures = 15;
    else
        numFeatures = 23;
    end
else
    error(['Unknown patch feature version: ' num2str(version)]);
end

[XSIZE, YSIZE, ZSIZE, junk] = size(dt_lowres);
pfmap = zeros(XSIZE, YSIZE, ZSIZE, numFeatures);
% The start and end indices just avoid potentially falling off the edge of
% the image when indexing the full neighbourhood.
parfor k=(n+1):(ZSIZE-n)
%for k=(n+1):(ZSIZE-n)
    for j=(n+1):(YSIZE-n)
        for i=(n+1):(XSIZE-n)
                
            ipatch = dt_lowres((i-n):(i+n),(j-n):(j+n),(k-n):(k+n),3:8);
            
            % Process if the central voxel is foreground
            if(dt_lowres(i,j,k,1)>=0)
    
                % Check whether any values in the patch are missing.
                isbg = (dt_lowres((i-n):(i+n),(j-n):(j+n),(k-n):(k+n),1)<0);
                if(sum(isbg(:))>0)
                    % Create a list of which elements of ipatch(:) are
                    % missing.
                    ismissing=repmat(isbg(:),[6,1]);
                    missing_inds = find(ismissing);
                    
                    % Construct submatrices of patch covariance.
                    T_11 = T(missing_inds, missing_inds);
                    T_22 = T(setdiff(1:length(ismissing), missing_inds), setdiff(1:length(ismissing), missing_inds));
                    T_12 = T(missing_inds, setdiff(1:length(ismissing), missing_inds));

                    % General expression for the conditional mean is (Prince page 48)
                    % \mu_1+\Sigma_{21}^T\Sigma_{22}^{-1}(x_2 - \mu_2). 
                    missing_cond_mean = comipatchmean(missing_inds)' + T_12*inv(T_22)*(ipatch(setdiff(1:length(ismissing), missing_inds))' - comipatchmean(setdiff(1:length(ismissing), missing_inds))');
                    
                    % Put the conditional mean values into the input patch.
                    ipatch(missing_inds) = missing_cond_mean;
                end
                
                if(version==1)
                    pfmap(i,j,k,:) = PatchFeatures(ipatch(:));
                elseif(version==2)
                    pfmap(i,j,k,:) = PatchFeaturesV2(ipatch(:));
                elseif(version==4)
                    pfmap(i,j,k,:) = PatchFeaturesV4(ipatch(:));
                elseif(version==5)
                    pfmap(i,j,k,:) = PatchFeaturesV5(ipatch(:));
                elseif(version==6)
                    pfmap(i,j,k,:) = PatchFeaturesV6(ipatch(:));
                else
                    error(['Unknown patch feature version: ' num2str(version)]);
                end
            end
        end
    end
end


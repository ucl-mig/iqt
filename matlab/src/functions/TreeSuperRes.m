function [dtrecon, accum, accum_var, leafmap, mdistmap] = TreeSuperRes(dt_lowres, tree, n, m, patch_feature_map, us, overlap, rescale_factor,  max_layer,  scale_const, input_patch_scaling, output_patch_scaling)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

if(nargin<7)
    overlap = 'ignore';
end
if(nargin<8)
    rescale_factor = 1.0;
end
if(nargin<9)
    max_layer = inf;
end
if(nargin<10)
    scale_const = 1E-3;
end
if(nargin<11)
    input_patch_scaling = 1;
end
if(nargin<12)
    output_patch_scaling = 1;
end

[XSIZE, YSIZE, ZSIZE, COMP] = size(dt_lowres);
dtrecon = zeros(XSIZE*us, YSIZE*us, ZSIZE*us, COMP);
dtrecon(:,:,:,1) = -1;
accum = zeros(XSIZE*us, YSIZE*us, ZSIZE*us);
accum_var = zeros(XSIZE*us, YSIZE*us, ZSIZE*us);
if(strcmp(overlap, 'cov_weighted_average'))
    accum = zeros(XSIZE*us, YSIZE*us, ZSIZE*us, COMP-2, COMP-2);
end
leafmap = zeros(XSIZE*us, YSIZE*us, ZSIZE*us);
mdistmap = zeros(XSIZE*us, YSIZE*us, ZSIZE*us);

% The start and end indices just avoid potentially falling off the edge of
% the image when indexing the full neighbourhood.
for k=(n+1):(ZSIZE-n)
    display(sprintf('Slice %i of %i.', k, ZSIZE));
    for j=(n+1):(YSIZE-n)
        for i=(n+1):(XSIZE-n)
                
            ipatch = dt_lowres((i-n):(i+n),(j-n):(j+n),(k-n):(k+n),3:COMP)*rescale_factor;
            
            % Process if the whole patch is foreground
            if(min(min(min(dt_lowres((i-n):(i+n),(j-n):(j+n),(k-n):(k+n),1))))>=0)
    
                % Find the right leaf
                features = patch_feature_map(i,j,k,:);
                tn = 1;
                count = 1; %track number of layers the input patch has traversed. 
                while(~tree{tn}.IsLeaf) && (count < max_layer)
                    if(features(tree{tn}.SplitFeatureIndex)<tree{tn}.FeatureThreshold)
                        tn = tree{tn}.LeftChildIndex;
                    else
                        tn = tree{tn}.RightChildIndex;
                    end
                    count = count + 1;
                end
                leafmap((us*(i-1)+1):(us*i),(us*(j-1)+1):(us*j),(us*(k-1)+1):(us*k)) = tn;
                
%                 % Compute Mahalanobis distance of input patch to training
%                 % patches for this leaf.
                 ipatch_scaled = (ipatch(:)')./input_patch_scaling;
%                 mean_diff = tree{tn}.PatchStats.ipatchmean - ipatch_scaled;
%                 % Mahalanobis distance
%                 % mdist = mean_diff*tree{tn}.PatchStats.iT*mean_diff';
%                 % Euclidean distance
%                 mdist = sqrt(sum(mean_diff.^2));
%                 mdistmap((us*(i-1)+1):(us*i),(us*(j-1)+1):(us*j),(us*(k-1)+1):(us*k)) = mdist;
                
                % Get corresponding output.
                opatchrecon = NodeRegress(tree{tn}, [ipatch_scaled scale_const]);%tree{tn}.Transformation.M'*ipatch(:);
                
                % Rescale it as specified
                opatchrecon = opatchrecon.*output_patch_scaling;
                
                % And recreate structure
                opatchrecon = reshape(opatchrecon, [m,m,m,COMP-2])/rescale_factor;

                if(strcmp(overlap, 'ignore'))
                    % Paste the central part of the reconstructed patch into
                    % the reconstructed image;
                    dtrecon((us*(i-1)+1):(us*i),(us*(j-1)+1):(us*j),(us*(k-1)+1):(us*k),3:COMP) = opatchrecon((1+(m-us)/2):((m+us)/2), (1+(m-us)/2):((m+us)/2), (1+(m-us)/2):((m+us)/2), :);

                    % Add in a best guess b=0 from original patch.  Here it
                    % just comes from the centre of the input patch. This could
                    % be done better...
                    logS0 = dt_lowres(i,j,k,2);
                    dtrecon((us*(i-1)+1):(us*i),(us*(j-1)+1):(us*j),(us*(k-1)+1):(us*k),2) = logS0;

                    % Label only reconstructed voxels as foreground.
                    dtrecon((us*(i-1)+1):(us*i),(us*(j-1)+1):(us*j),(us*(k-1)+1):(us*k),1) = 0;
                    
                elseif(strcmp(overlap, 'average'))
                    % Add patch to the reconstructed image;
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) + opatchrecon;

                    % Add best guess b=0 from original patch.
                    logS0 = dt_lowres(i,j,k,2);
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) + logS0;

                    % Increment the accumulator to keep track of
                    % normalization required for average.
                    accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2)) = accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2))+1;
                    
                    % Label only reconstructed voxels as foreground.
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),1) = 0;
                elseif(strcmp(overlap, 'weighted_average_valid'))
                    % Add patch to the reconstructed image weighted by
                    % local variance estimate
                    local_var = mean(diag(tree{tn}.Transformation.Sv));
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) + opatchrecon/local_var;

                    % Add best guess b=0 from original patch.
                    logS0 = dt_lowres(i,j,k,2);
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) + logS0;

                    % Add one over the local variance estimate to the accumulator to keep track of
                    % normalization required for average.
                    accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2)) = accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2))+1/local_var;
                    accum_var((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2)) = accum_var((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2))+ local_var;

                    % Label only reconstructed voxels as foreground.
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),1) = 0;
                elseif(strcmp(overlap, 'weighted_average'))
                    % Add patch to the reconstructed image weighted by
                    % local variance estimate
                    local_var = mean(diag(tree{tn}.Transformation.S));
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) + opatchrecon/local_var;

                    % Add best guess b=0 from original patch.
                    logS0 = dt_lowres(i,j,k,2);
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) + logS0;

                    % Add one over the local variance estimate to the accumulator to keep track of
                    % normalization required for average.
                    accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2)) = accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2))+1/local_var;
                    accum_var((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2)) = accum_var((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2))+ local_var;

                    % Label only reconstructed voxels as foreground.
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),1) = 0;
                
                elseif(strcmp(overlap, 'weighted_average_iso'))
                    % Add patch to the reconstructed image weighted by
                    % local variance estimate on validation data for
                    % isotropic trees.
                    local_var = tree{tn}.Transformation.Sv;
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) + opatchrecon/local_var;

                    % Add best guess b=0 from original patch.
                    logS0 = dt_lowres(i,j,k,2);
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) + logS0;

                    % Add one over the local variance estimate to the accumulator to keep track of
                    % normalization required for average.
                    accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2)) = accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2))+1/local_var;
                    accum_var((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2)) = accum_var((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2))+ local_var;

                    % Label only reconstructed voxels as foreground.
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),1) = 0;

                elseif(strcmp(overlap, 'cov_weighted_average'))
                    
                    % Calculate the contribution to each subvoxel.
                    covind = 1;
                    for kp=1:m
                        for jp=1:m
                            for ip=1:m
                                icov = inv(tree{tn}.Transformation.Sv(covind:(m^3):end, covind:(m^3):end));
                                dtcontrib = icov*opatchrecon(covind:(m^3):end)';
                                
                                % Add it to the recon image
                                dtrecon(us*(i-1)-(m-us)/2 + ip, us*(j-1)-(m-us)/2 + jp, us*(k-1)-(m-us)/2 + kp, 3:COMP) = squeeze(dtrecon(us*(i-1)-(m-us)/2 + ip, us*(j-1)-(m-us)/2 + jp, us*(k-1)-(m-us)/2 + kp, 3:COMP)) + dtcontrib;
                                
                                % Increment the normalisation accumulator.
                                accum(us*(i-1)-(m-us)/2 + ip, us*(j-1)-(m-us)/2 + jp, us*(k-1)-(m-us)/2 + kp, :,:) = squeeze(accum(us*(i-1)-(m-us)/2 + ip, us*(j-1)-(m-us)/2 + jp, us*(k-1)-(m-us)/2 + kp, :,:)) + icov;
                                
                                % Increment the patch element index.
                                covind = covind + 1;
                            end
                        end
                    end
                    
                    % Add best guess b=0 from original patch.
                    logS0 = dt_lowres(i,j,k,2);
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) = logS0;
                    
                    % Label only reconstructed voxels as foreground.
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),1) = 0;
                else
                    error(sprintf('Unknown overlap strategy: %s.', overlap));
                end
            end
        end
    end
end

% Renormalize
if(strcmp(overlap, 'average') || strcmp(overlap, 'weighted_average') || strcmp(overlap, 'weighted_average_valid')||strcmp(overlap, 'weighted_average_iso'))
    fg_inds = find(dtrecon(:,:,:,1)==0);
    for i=2:COMP
        v = dtrecon(:,:,:,i);
        v(fg_inds) = v(fg_inds)./accum(fg_inds);
        dtrecon(:,:,:,i) = v;
    end
end




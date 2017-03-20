function [dtrecon, accum, leafmap, mx_hash] = rt_TreeSuperResEdge(dt_lowres, tree, n, m, patch_feature_map, us, T, comipatchmean, overlap, rescale_factor, scale_const, input_patch_scaling, output_patch_scaling, mx_hash)

if(nargin<9)
    overlap = 'ignore';
end
if(nargin<10)
    rescale_factor = 1.0;
end
if(nargin<11)
    scale_const = 1E-3;
end
if(nargin<12)
    input_patch_scaling = 1;
end
if(nargin<13)
    output_patch_scaling = 1;
end
if(nargin<14)
    mx_hash.n = n;
    mx_hash.m = m;
    mx_hash.T = T;
    mx_hash.comipatchmean = comipatchmean;
end

% This is a good limit using double precision matrices on abner-7-1, which
% has 128G RAM.
%max_hash_size = 20000;
% Changing to storing matrices with single precision, we can double the
% size, which then covers whole image with n=2.
max_hash_size = 40000;

[XSIZE, YSIZE, ZSIZE, COMP] = size(dt_lowres);
dtrecon = zeros(XSIZE*us, YSIZE*us, ZSIZE*us, COMP);
dtrecon(:,:,:,1) = -1;
accum = zeros(XSIZE*us, YSIZE*us, ZSIZE*us);
if(strcmp(overlap, 'cov_weighted_average'))
    accum = zeros(XSIZE*us, YSIZE*us, ZSIZE*us, COMP-2, COMP-2);
end
leafmap = zeros(XSIZE*us, YSIZE*us, ZSIZE*us);

% The start and end indices just avoid potentially falling off the edge of
% the image when indexing the full neighbourhood.
repeats = 0;
warning_given = 0;
for k=(n+1):(ZSIZE-n)
    display([sprintf('Slice %i of %i. ', k, ZSIZE) 'Edge stats: Repeats: ' num2str(repeats) '; Stored fields: ' num2str(length(fields(mx_hash)))]);
    for j=(n+1):(YSIZE-n)
        for i=(n+1):(XSIZE-n)
                
            ipatch = dt_lowres((i-n):(i+n),(j-n):(j+n),(k-n):(k+n),3:COMP);
            
            % Process if the central voxel is foreground
            if(min(min(min(dt_lowres(i,j,k,1))))>=0)
    
                % Check whether any values in the patch are missing.
                isbg = (dt_lowres((i-n):(i+n),(j-n):(j+n),(k-n):(k+n),1)<0);
                if(sum(isbg(:))>0)
                    
                    % Construct unique identifier for this pattern of
                    % missing values. It converts the binary number with
                    % all consecutive entries of isbg into hex. Note that
                    % bin2dec will process a max number length of 52 so we 
                    % process in chunks of 48.
                    s = '';
                    for bl = 1:length(isbg(:))
                        s = [s num2str(isbg(bl))];
                    end
                    bgid = 'Z';
                    for bl = 1:(fix(length(s)/48)+1)
                        bgid = [bgid dec2hex(bin2dec(s(((bl-1)*48+1):min([length(s), bl*48]))))];
                    end
                    
                    % Create a list of which elements of ipatch(:) are
                    % missing.
                    ismissing=repmat(isbg(:),[COMP-2,1]);
                    missing_inds = find(ismissing);
                    
                    % Compute required matrix and store in hash table.
                    if(isfield(mx_hash, bgid))
                        T_12divT_22 = eval(['mx_hash.' bgid]);
                        repeats = repeats + 1;
                    else                        
                        % Construct submatrices of patch covariance.
                        %T_11 = T(missing_inds, missing_inds);
                        T_22 = T(setdiff(1:length(ismissing), missing_inds), setdiff(1:length(ismissing), missing_inds));
                        T_12 = T(missing_inds, setdiff(1:length(ismissing), missing_inds));

                        T_12divT_22 = single(T_12*pinv(T_22));
                        if(length(fields(mx_hash))<max_hash_size)
                            eval(['mx_hash.' bgid ' = T_12divT_22;']);
                        else
                            if(~warning_given)
                                warning('Reached max matrix hash size');
                                warning_given = 1;
                            end
                        end
                    end
                    
                    % General expression for the conditional mean is (Prince page 48)
                    % \mu_1+\Sigma_{21}^T\Sigma_{22}^{-1}(x_2 - \mu_2). 
                    missing_cond_mean = comipatchmean(missing_inds)' + T_12divT_22*(ipatch(setdiff(1:length(ismissing), missing_inds))' - comipatchmean(setdiff(1:length(ismissing), missing_inds))');
                    
                    % Put the conditional mean values into the input patch.
                    ipatch(missing_inds) = missing_cond_mean;
                end
                
                % Rescale
                ipatch = ipatch*rescale_factor;
                
                features = patch_feature_map(i,j,k,:);
                tn = 1;
                while(~tree{tn}.IsLeaf)
                    if(features(tree{tn}.SplitFeatureIndex)<tree{tn}.FeatureThreshold)
                        tn = tree{tn}.LeftChildIndex;
                    else
                        tn = tree{tn}.RightChildIndex;
                    end
                end
                leafmap(i,j,k) = tn;
                
                % vectorise and add a constant:
                ipatch_scaled = (ipatch(:)')./input_patch_scaling; %(!): in this c
                ipatch_scaled = [ipatch_scaled, scale_const];
                
                if (strcmp(overlap, 'std_bay_weighted_average'))
                    ipatch_scaled = (ipatch_scaled - tree{tn}.Transformation.train_mean)./tree{tn}.Transformation.train_std;
                    opatchrecon = rt_BayNodeRegress(tree{tn},ipatch_scaled',[],[],'MAP');%tree{tn}.Transformation.M'*ipatch(:);
                    %opatchrecon = NodeRegress(tree{tn}, [ipatch_scaled scale_const]);%tree{tn}.Transformation.M'*ipatch(:);
                else
                    % Use linear mapping to get corresponding output.
                    opatchrecon = NodeRegress(tree{tn}, ipatch_scaled);%tree{tn}.Transformation.M'*ipatch(:);
                end 
                
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
                elseif(strcmp(overlap, 'weighted_average'))
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
                    
                    % Label only reconstructed voxels as foreground.
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),1) = 0;
                
                 elseif (strcmp(overlap, 'std_bay_weighted_average'))
                    % Add patch to the reconstructed image weighted by
                    % local variance estimate
                    local_var = ipatch_scaled*tree{tn}.Transformation.A*ipatch_scaled' + tree{tn}.HyperPara.var_noise;
                    %local_var = mean(diag(tree{tn}.Transformation.Sv));
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),3:COMP) + opatchrecon/local_var;

                    % Add best guess b=0 from original patch (?) - why is this necessary?
                    logS0 = dt_lowres(i,j,k,2);
                    dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) = dtrecon((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2),2) + logS0;

                    % Add one over the local variance estimate to the accumulator to keep track of
                    % normalization required for average.
                    accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2)) = accum((us*(i-1)+1-(m-us)/2):(us*i+(m-us)/2),(us*(j-1)+1-(m-us)/2):(us*j+(m-us)/2),(us*(k-1)+1-(m-us)/2):(us*k+(m-us)/2))+1/local_var;

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
if(strcmp(overlap, 'average') || strcmp(overlap, 'weighted_average') || strcmp(overlap, 'std_bay_weighted_average'))
    fg_inds = find(dtrecon(:,:,:,1)==0);
    for i=2:COMP
        v = dtrecon(:,:,:,i);
        v(fg_inds) = v(fg_inds)./accum(fg_inds);
        dtrecon(:,:,:,i) = v;
    end
end




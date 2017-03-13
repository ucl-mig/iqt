function exec_BigbassCreateTrainingSetPatchLibs(fi, OutputPath, data_folders, SubPath, dt_name)

    output_dir = [OutputPath '/'  data_folders{fi} '/' SubPath ];
    
    disp( ' ' );
    disp( [ 'PROCESSING : ' data_folders{fi} ] );
    disp( ' ' );

    % The downsampling factor for the input.
    for ds = 2:4;

        disp(sprintf('Downsampling rate of : %i', ds));

        %create_feature_vectors( path_data );

        dt_lowres_name = sprintf('%slowres_%i_', dt_name, ds) ;

        dt_suffix = '' ;
        %dt_suffix = '_log' ;

        % load the DT images

        [mask,elsp] = read_std_nii([ output_dir dt_name '1.nii' ]);

        dim1 = size(mask,1);
        dim2 = size(mask,2);
        dim3 = size(mask,3);

        % The size of the subvoxel array (output) is mxmxm
        m = ds;
        % radius of low-res neighbourhood, ie input is radius_lr_patchxradius_lr_patchxradius_lr_patch
        for radius_lr_patch = 1:3

            tic();

            inds = find(mask==0);

            disp(sprintf('computing masks for n=%i.', radius_lr_patch));

            %--------------------------------------------------------------------------
            % compute area in which the MxMxM-blocks are computed only on valid voxels
            % the valid points are marked by setting the first voxel of the MxMxM block to 0/1 

            mask_lr_valid_blocks = mask*0;

            for x=inds'
            % for k=1:dim3
            %     for j=1:dim2
            %         for i=1:dim1

                        [i,j,k] = ind2sub( [dim1,dim2,dim3] , x );

                        if ( i+m-1<=dim1 && j+m-1<=dim2 && k+m-1<=dim3 )

                            % if one singe voxel from the MxMxM-block is invalid, do not use the block. 
                            tmp = mask( i:i+m-1 , j:j+m-1 , k:k+m-1 );
                            mask_lr_valid_blocks( i,j,k ) = min( tmp(:)==0 );

                        end

            %         end
            %     end
            % end
            end



            %--------------------------------------------------------------------------
            % compute blocks which are valid for reconstruction.
            % For this, all blocks of a patch around a valid MxMxM-block to be reconstructed
            % must also be valid since they are used for the estimation.
            % the valid points are marked by setting the first voxel of the MxMxM block of interest (ie. block to be reconstructed) to 0/1 
            mask_hr_valid_blocks_for_rec = mask*0; % valid only if surrounding is also valid

            for x=inds'
            % for k=1:dim3
            %     for j=1:dim2
            %         for i=1:dim1

                        [i,j,k] = ind2sub( [dim1,dim2,dim3] , x );

                        if ( i+(radius_lr_patch+1)*ds-1<=dim1 && j+(radius_lr_patch+1)*ds-1<=dim2 && k+(radius_lr_patch+1)*ds-1<=dim3 && ...
                             i-radius_lr_patch*ds>=1 && j-radius_lr_patch*ds>=1 && k-radius_lr_patch*ds>=1 )

                            tmp = mask( i-radius_lr_patch*ds:i+(radius_lr_patch+1)*ds-1 , j-radius_lr_patch*ds:j+(radius_lr_patch+1)*ds-1 , k-radius_lr_patch*ds:k+(radius_lr_patch+1)*ds-1 );
                            mask_hr_valid_blocks_for_rec( i,j,k ) = min( tmp(:)==0 );

                        end

            %         end
            %     end
            % end
            end

            indices_valid_features = find(mask_hr_valid_blocks_for_rec>0);


            %--------------------------------------------------------------------------
            % collect valid feature vectors

            if strcmp(dt_suffix,'')
                indices_patchlib_to_volume = indices_valid_features;
                save( [ output_dir sprintf('indices_patchlib_to_volumeDS%02i_N%02i', ds, radius_lr_patch) ] , 'indices_patchlib_to_volume' );
            end


            % Low-Res -------------------------------

            disp('creating patch library LOW-RES...')

            dt_lr = zeros( dim1 , dim2 , dim3 , 6 );
            for i=1:6

                dt_lr(:,:,:,i) = read_std_nii([ output_dir dt_lowres_name num2str(i+2) dt_suffix '.nii' ]);

            end

            ipatchlib = zeros( numel(indices_valid_features) , 6*(2*radius_lr_patch+1)^3, 'single' );

            for ind_i=1:numel(indices_valid_features)

                ind = indices_valid_features(ind_i);

                [i,j,k] = ind2sub( [dim1,dim2,dim3] , ind );

                ipatch = dt_lr( i-radius_lr_patch*ds:ds:i+(radius_lr_patch+1)*ds-1 , j-radius_lr_patch*ds:ds:j+(radius_lr_patch+1)*ds-1 , k-radius_lr_patch*ds:ds:k+(radius_lr_patch+1)*ds-1 , : );

                ipatchlib(ind_i,:) = ipatch(:);

            end
            disp('done.')

            disp('saving patch library...')
            save( [ output_dir sprintf('ipatchlibDS%02i_N%02i', ds, radius_lr_patch) dt_suffix ] , 'ipatchlib', '-v7.3');
            disp('done.')

%             clear ipatchlib
%             clear dt_lr
            ipatchlib = [];
            dt_lr = [];


            if(radius_lr_patch==1)

                % High-Res ------------------------------

                disp('creating patch library HIGH-RES...')

                dt = zeros( dim1 , dim2 , dim3 , 6 );
                for i=1:6

                    dt(:,:,:,i) = read_std_nii([ output_dir dt_name num2str(i+2) dt_suffix '.nii' ]);

                end

                % Create three output patch libraries.  The first includes just the
                % central block of voxels, ie dsxdsxds patch.  The second includes a
                % shell of voxels outside that (ds+2)x(ds+2)x(ds+2).  The last two
                % shells: (ds+4)x(ds+4)x(ds+4).
                opatchlib_ds = zeros( numel(indices_valid_features) , 6*ds^3 , 'single');
                if(ds==2 && radius_lr_patch==1)
                    opatchlib_dsp1 = zeros( numel(indices_valid_features) , 6*(ds+2)^3 , 'single');
                    opatchlib_dsp2 = zeros( numel(indices_valid_features) , 6*(ds+4)^3, 'single' );
                end

                for ind_i=1:numel(indices_valid_features)

                    ind = indices_valid_features(ind_i);

                    [i,j,k] = ind2sub( [dim1,dim2,dim3] , ind );

                    opatch = dt( i:i+ds-1 , j:j+ds-1 , k:k+ds-1 , : );
                    opatchlib_ds(ind_i,:) = opatch(:);

                    if(ds==2 && radius_lr_patch==1)
                        opatch = dt( (i-1):(i+ds) , (j-1):(j+ds) , (k-1):(k+ds) , : );
                        opatchlib_dsp1(ind_i,:) = opatch(:);

                        opatch = dt( (i-2):(i+ds+1) , (j-2):(j+ds+1) , (k-2):(k+ds+1) , : );
                        opatchlib_dsp2(ind_i,:) = opatch(:);
                    end

                end
                disp('done.')

                disp('saving patch library...')
                opatchlib = opatchlib_ds;
                save( [ output_dir sprintf('opatchlibDS%02i_N%02i_M%02i', ds, radius_lr_patch, ds) dt_suffix ] , 'opatchlib', '-v7.3' );
                if(ds==2 && radius_lr_patch==1)
                    opatchlib = opatchlib_dsp1;
                    save( [ output_dir sprintf('opatchlibDS%02i_N%02i_M%02i', ds, radius_lr_patch, ds+2) dt_suffix ] , 'opatchlib', '-v7.3' );
                    opatchlib = opatchlib_dsp2;
                    save( [ output_dir sprintf('opatchlibDS%02i_N%02i_M%02i', ds, radius_lr_patch, ds+4) dt_suffix ] , 'opatchlib', '-v7.3' );
                end
                disp('done.')

%                 clear opatchlib opatchlib_ds opatchlib_dsp1 opatchlib_dsp2
%                 clear dt
                opatchlib = [];
                opatchlib_ds = [];
                opatchlib_dsp1 =[];
                opatchlib_dsp2 = [];
                dt = [];

                indices_patchlib_to_volumeN1 = indices_patchlib_to_volume;

            else
                indices_to_N1_indices = zeros(length(indices_patchlib_to_volume),1);
                ind = 1;
                for i=1:length(indices_patchlib_to_volume)
                    while(indices_patchlib_to_volumeN1(ind) ~= indices_patchlib_to_volume(i))
                        ind = ind + 1;
                    end
                    indices_to_N1_indices(i) = ind;
                end

                save( [ output_dir sprintf('indices_to_N1_indicesDS%02i_N%02i', ds, radius_lr_patch) dt_suffix ] , 'indices_to_N1_indices', '-v7.3' );
            end


            disp([ 'duration: took ' num2str(toc) ' seconds.' ]);
        end
    end


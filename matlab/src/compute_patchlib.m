function compute_patchlib(input_dir, output_dir, data_folders, params)

% Fetch the parameters:
dt_name = params.dt_name;
ds = params.upsample_rate; % downsampling rate
input_radius = params.input_radius; % the radius of the low-res patch. 
output_radius = ds;

% compute patch libaries sequentiatlly for all subjects.
for fi = 1:length(data_folders)
    if(~exist([output_dir '/'  data_folders{fi} ]))
        mkdir([output_dir '/'  data_folders{fi} ]);
    end
   
    input_folder = [input_dir '/' data_folders{fi}];
    output_folder = [output_dir '/' data_folders{fi}];
    
    disp( ' ' );
    disp( [ 'Creating a patch library from: ' data_folders{fi} ] );
    disp( ' ' );

    % The downsampling factor for the input.
    disp(sprintf('Downsampling rate of : %i', ds));
    dt_lowres_name = sprintf('%slowres_%i_', dt_name, ds);
    
    % radius of low-res neighbourhood;
    
    %----------- Find the indices of valid patches ------------------------
    % We only use patches fully contained inside the brain in both high-res
    % and low-res space. We refer to such patches as valida patches.
    % Here the location of the central voxel of valida patches are computed
    % in the high-res space.
    % background (-1) or brain (0) voxels in the mask.
    
    disp(sprintf(['Computing locations of valid patch pairs: '...
                  'input(%ix%ix%i)  =>  output(%ix%ix%i)'],...
        2*input_radius+1,2*input_radius+1,2*input_radius+1,...
        ds,ds,ds));
    
    % load the DT brain mask: 
    [mask, elsp] = read_std_nii([ input_folder '/' dt_name '1.nii' ]);
    dim1 = size(mask,1); dim2 = size(mask,2); dim3 = size(mask,3);
    
    mask_hr_valid_blocks_for_rec = mask*0;
    inds = find(mask==0);
    
    tic()
    for x=inds'
        
        [i,j,k] = ind2sub( [dim1,dim2,dim3] , x );
        
        if ( i+(input_radius+1)*ds-1<=dim1 && ...
             j+(input_radius+1)*ds-1<=dim2 && ...
             k+(input_radius+1)*ds-1<=dim3 && ...
             i-input_radius*ds>=1 && j-input_radius*ds>=1 && ...
             k-input_radius*ds>=1 )

            tmp = mask( i-input_radius*ds:i+(input_radius+1)*ds-1,...
                        j-input_radius*ds:j+(input_radius+1)*ds-1,...
                        k-input_radius*ds:k+(input_radius+1)*ds-1 );
                    
            mask_hr_valid_blocks_for_rec( i,j,k ) = min( tmp(:)==0 );
        end
    end

    % collect the indices of valid feature vectors:
    indices_valid_features = find(mask_hr_valid_blocks_for_rec>0);
    
    % Save the indices (not neccesary).
    indices_patchlib_to_volume = indices_valid_features;
    save( [ output_folder sprintf('indices_patchlib_to_volumeDS%02i_N%02i',...
            ds, input_radius) ] , 'indices_patchlib_to_volume' );

    % ------------ Extract patches from low-res DTI (input) --------------
    % Extract and rasterise all valid patches from low-res DTI, and put
    % them all in a big single matrix. 
    % Each patch is of size (2*input_radius+1)^3*6
    % Save as a mat file.
    
    disp('creating patch library LOW-RES...')
    dt_lr = zeros( dim1 , dim2 , dim3 , 6 );
    for i=1:6
        dt_lr(:,:,:,i) = read_std_nii(...
                         [input_folder ...
                          dt_lowres_name num2str(i+2) '.nii' ]);
    end
    
    % define the input patch library.
    ipatchlib = zeros(numel(indices_valid_features), 6*(2*input_radius+1)^3,...
                     'single' );

    for ind_i=1:numel(indices_valid_features)

        ind = indices_valid_features(ind_i);

        [i,j,k] = ind2sub( [dim1,dim2,dim3] , ind );

        ipatch = dt_lr( i-input_radius*ds:ds:i+(input_radius+1)*ds-1,...
                        j-input_radius*ds:ds:j+(input_radius+1)*ds-1,...
                        k-input_radius*ds:ds:k+(input_radius+1)*ds-1 , : );

        ipatchlib(ind_i,:) = ipatch(:);

    end
    disp('done.')

    disp('saving patch library...')
    save([output_folder ...
          sprintf('ipatchlibDS%02i_N%02i', ds, input_radius)], ...
          'ipatchlib', '-v7.3');
    disp('done.')

    ipatchlib = [];
    dt_lr = [];
    
    % ------------ Extract patches from high-res DTI (output) -------------    
    % Each high-res patch is of size ds^3*6
    disp('creating patch library HIGH-RES...')
    dt = zeros( dim1 , dim2 , dim3 , 6 );
    
    % load the high-res dti
    for i=1:6
        dt(:,:,:,i) = read_std_nii([input_folder ...
                                    dt_name num2str(i+2), '.nii' ]);
    end
    
    % Define the input patch library.
    opatchlib = zeros( numel(indices_valid_features) , 6*ds^3 , 'single');

    for ind_i=1:numel(indices_valid_features)

        ind = indices_valid_features(ind_i);

        [i,j,k] = ind2sub( [dim1,dim2,dim3] , ind );

        opatch = dt( i:i+ds-1 , j:j+ds-1 , k:k+ds-1 , : );
        opatchlib(ind_i,:) = opatch(:);
    end
    disp('done.')

    disp('saving patch library...')
    save([output_folder ...
         sprintf('opatchlibDS%02i_N%02i_M%02i', ds, input_radius, ds)],...
         'opatchlib', '-v7.3' );
    disp('done.')

    % clear for memory.
    opatchlib = [];
    dt = [];

    disp([ 'duration: took ' num2str(toc) ' seconds.' ]);
end
    

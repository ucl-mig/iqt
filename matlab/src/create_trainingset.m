function create_trainingset(input_dir, output_dir, data_folders, params)
% Merge patch libraries from a training set of data sets.

SampleRate = params.subsample_rate;
rnds = params.no_rnds;
n = params.input_radius; % the radius of the low-res patch. 
ds = params.upsample_rate;
m = ds;

% Loop over randomizations
for si = 1:rnds
    for i=1:length(data_folders)
        load([input_dir '/' data_folders{i} '/' ...
              sprintf('ipatchlibDS%02i_N%02i', ds, n)]);
        load([input_dir '/' data_folders{i} '/' ...
              sprintf('opatchlibDS%02i_N01_M%02i', ds, m)]);
          
%         if(n>1)
%             % Need to subsample the opatchlib to match the ipatchlib
%             load([intput_dir '/' data_folders{i} sprintf('indices_to_N1_indicesDS%02i_N%02i', ds, n)]);
%             opatchlib = opatchlib(indices_to_N1_indices,:);
%         end

        keep = find(rand(1,length(ipatchlib))<(1/SampleRate));
        if(i==1)
            comipatchlib = ipatchlib(keep,:);
            comopatchlib = opatchlib(keep,:);
        else
            comipatchlib = [comipatchlib; ipatchlib(keep,:)];
            comopatchlib = [comopatchlib; opatchlib(keep,:)];
        end
        
        % Keep track of which voxels are included.
        patchlibindices{i} = keep;
    end

    % Save some space...
    comipatchlib = single(comipatchlib);
    comopatchlib = single(comopatchlib);

    filename = sprintf('PatchLibs_DS%02i_%ix%ix%i_%ix%ix%i_Subi%03i_%04i.mat', ds, 2*n+1,2*n+1,2*n+1, m, m, m, SampleRate, si);
    disp(['Saving the patch-library as' filename ])
    save([output_dir '/' filename], 'comipatchlib', 'comopatchlib', 'patchlibindices', '-v7.3');
end



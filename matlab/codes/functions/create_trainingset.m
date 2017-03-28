%% Create training sets from the patch-libraries of selected subjects.
% We randomly subsample from respective patch-libraries and merge them to
% form a training set. 

function create_trainingset(input_dir, output_dir, data_folders, settings)


SampleRate = settings.subsample_rate;
rnds = settings.no_rnds;
n = settings.input_radius; % the radius of the low-res patch. 
ds = settings.upsample_rate;
m = ds;

if(~exist(output_dir))
    mkdir(output_dir);
end

% Loop over randomizations
for si = 1:rnds
    for i=1:length(data_folders)
        load([input_dir '/' data_folders{i} '/' ...
              sprintf('ipatchlibDS%02i_N%02i', ds, n)]);
        load([input_dir '/' data_folders{i} '/' ...
              sprintf('opatchlibDS%02i_N%02i_M%02i', ds, n, m)]);
          
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

    filename = sprintf('PatchLibs_DS%02i_%ix%ix%i_%ix%ix%i_Sub%03i_%04i.mat', ds, 2*n+1,2*n+1,2*n+1, m, m, m, SampleRate, si);
    disp(['Saving the patch-library as: ' filename ])
    save([output_dir '/' filename], 'comipatchlib', 'comopatchlib', 'patchlibindices', '-v7.3');
end
disp(['The training set is available at: ' output_dir])



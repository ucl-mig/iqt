function visualise_results(hires_dir, lowres_dir, recon_dir, ...
                           data_folders, sub_path, dt_pref, ...
                           ds_rate, no_rnds, sample_rate, ...
                           input_radius, fv, edge_recon, ...
                           flip_dim)
% RECONSTRUCT_RANDOMFORESTS
%   Computes IQT super-resolution DTI from the input DTI (low-resolution).
%
%   Args: 
%       HIRES_DIR: Input root dir of hi-res DTI (='' if not available)
%       LOWRES_DIR: Output root dir of low-res DTI
%       RECON_DIR: Root dir containing the RF reconstructed DTI
%       DATA_FOLDERS: Subjects to process (eg HCP subjects)
%       SUB_PATH: Input/Output data access is done as:
%                  input_dir/data_folder/sub_path/file_name
%                  output_dir/data_folder/sub_path/file_name
%       DT_PREF: DTI output prefix
%       DS_RATE: Super-resolution factor
%       NO_RNDS: number of trees in the RF
%       SAMPLE_RATE: data sampling rate for random sub-sampling
%       INPUT_RADIUS: the input is a cubic patch of size (2*INPUT_RADIUS+1)^3
%       FV: feature version used for computing features for tree training
%       EDGE_RECON: flag that indicates whether to reconstructs edge or not
%                   (edge reconstruction is slow!)
%       FLIP_DIM (optional): The internal IQT code flips all images along 
%                            axis-1, due to historic reasons. But this 
%                            should not be the case for user data. 
%                            Normally, this flag should be off (default 0).
%
%   (always end directory paths with a forward/back slash)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

check_path(hires_dir);
check_path(lowres_dir);
check_path(recon_dir);
check_path(sub_path);

rescale_factor = 1.0;
scale_const = 1E-3;

ds=ds_rate;
n=input_radius;
m=ds;
trees_list = 1:no_rnds;
tail_name = sprintf('DS%02i_%ix%ix%i_%ix%ix%i_Sub%03i', ds, 2*n+1,2*n+1,2*n+1, m, m, m, sample_rate);

if edge_recon % Reconstruct on the boundary.
   % reconstruct and save the estimated DTI and precision map:
   output_subdir  = sprintf(['RF_Edge_V' int2str(fv) '_NoTree%02i_' tail_name '/'], length(trees_list));
else
   % reconstruct and save the estimated DTI and precision map:
   output_subdir  = sprintf(['RF_V' int2str(fv) '_NoTree%02i_' tail_name '/'], length(trees_list));
end

if ~exist('flip_dim', 'var')
    flip_dim = 0;
end

% Load low-res, high-res, estimate
for dataid = 1:length(data_folders)
    fprintf(['\nVisualising: ' data_folders{dataid} '\n']);
    hires_folder  = [hires_dir data_folders{dataid} '/' sub_path];
    lowres_folder = [lowres_dir data_folders{dataid} '/' sub_path];
    recon_folder  = [recon_dir data_folders{dataid} '/' sub_path output_subdir];
    
    % Load the data:
    fprintf('Loading low-res DTI...\n');
    file_low = [ lowres_folder dt_pref sprintf('lowres_%i_', ds)];
    dt_lr = ReadDT_Volume(file_low);
    if flip_dim>0
        dt_lr = dt_lr(1:ds:end,1:ds:end,1:ds:end,:);
    end
    
    fprintf('Loading reconstructed DTI...\n');
    file_est = [recon_folder dt_pref 'recon_'];
    dt_est = ReadDT_Volume(file_est);
    
    file_orig = [hires_folder dt_pref];
    if strcmp(hires_dir, '')
        dt_hr = [];
    else
        fprintf('Loading hi-res DTI...\n');
        dt_hr = ReadDT_Volume(file_orig);
    end
        
    % Take a slice, and compute MD, FA, CFA.
    zN = round(size(dt_est, 3)/2);
    fprintf('Taking z-slice %i\n', zN);
    slice_est = dt_est(:,:,zN,:); dt_est = [];
    slice_lr  = dt_lr(:,:,round(zN/ds),:); dt_lr = [];
    if ~isempty(dt_hr)
        slice_hr= dt_hr(:,:,zN,:); dt_hr = [];
    else
        slice_hr = [];
    end
    
    fprintf('Computing MD/FA/CFA ...\n');
    [md_lr, fa_lr, cfa_lr] = compute_MD_FA_CFA(slice_lr);
    [md_est, fa_est, cfa_est] = compute_MD_FA_CFA(slice_est);
    if ~isempty(slice_hr)
        [md_hr, fa_hr, cfa_hr] = compute_MD_FA_CFA(slice_hr);
    else
        md_hr = single(zeros(size(md_est)));
        fa_hr = single(zeros(size(fa_est)));
        cfa_hr = single(zeros(size(cfa_est)));
    end
    
    % Flip:
    md_lr = flipud(md_lr');
    md_est = flipud(md_est');
    md_hr = flipud(md_hr');
    
    fa_lr = flipud(fa_lr');
    fa_est = flipud(fa_est');
    fa_hr = flipud(fa_hr');
    
    cfa_lr2(:,:,1) = flipud(cfa_lr(:,:,1)');
    cfa_lr2(:,:,2) = flipud(cfa_lr(:,:,2)');
    cfa_lr2(:,:,3) = flipud(cfa_lr(:,:,3)');
    cfa_est2(:,:,1) = flipud(cfa_est(:,:,1)');
    cfa_est2(:,:,2) = flipud(cfa_est(:,:,2)');
    cfa_est2(:,:,3) = flipud(cfa_est(:,:,3)');
    cfa_hr2(:,:,1) = flipud(cfa_hr(:,:,1)');
    cfa_hr2(:,:,2) = flipud(cfa_hr(:,:,2)');
    cfa_hr2(:,:,3) = flipud(cfa_hr(:,:,3)');
    
    % Plot:
    margin = [0.02,0.02];
    fig=figure; 
    subplot_tight(3,3,1,margin)
    imshow(md_lr,[]);
    title('Low-res input')
    ylabel('MD')
    subplot_tight(3,3,2,margin)
    imshow(md_est,[]);
    title('IQT-RF outut')
    subplot_tight(3,3,3,margin)
    imshow(md_hr,[]);
    title('Ground truth high-res')
    
    subplot_tight(3,3,4,margin)
    imshow(cfa_lr2,[]);
    ylabel('CFA')
    subplot_tight(3,3,5,margin)
    imshow(cfa_est2,[]);
    subplot_tight(3,3,6,margin)
    imshow(cfa_hr2,[]);
    
    subplot_tight(3,3,7,margin)
    imshow(fa_lr,[]);
    ylabel('FA')
    subplot_tight(3,3,8,margin)
    imshow(fa_est,[]);
    subplot_tight(3,3,9,margin)
    imshow(fa_hr,[]);
    
    %Save as a .fig file.
    disp('Save the figure as a FIG file:')
    filename = [recon_folder 'image.fig'];
    disp(['see ' filename])
    saveas(fig,filename)
end

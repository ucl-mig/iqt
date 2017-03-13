function BigbassPrepareTrainingSet(SampleRate, n, tsid, ts, rnds, ds, m)
% Merge patch libraries from a training set of data sets.

if(nargin<2)
    n = 1;
end

if(nargin<3)
    tsid = 'Orig';
end
ExtraPath = '/T1w/Diffusion/';
if(strcmp(tsid, 'Orig'))
    % Original 8 random training set elements
    TrainingSet = {'101915','106016','120111','122317','130316','148335','153025','159340', '351938',  '390645',  '545345',  '586460',  '705341',  '749361',  '765056',  '951457'};
elseif(strcmp(tsid, 'Focus'))
    % 8 22-27 white non-hispanic right-handed females training set.
    TrainingSet = {'153025',  '103414',  '727654',  '151223',  '167743', '100307', '130316', '113215'};
elseif(strcmp(tsid, 'Diverse'))
    % Diverse training set.
    TrainingSet = {'992774', '125525', '205119', '133928', '570243', '448347', '654754', '153025'};
elseif(strcmp(tsid, 'Random1'))
    % Random training set.
    TrainingSet = {'106319', '117122', '133827', '140824', '158540', '196750', '205826', '366446'};
elseif(strcmp(tsid, 'Random2'))
    % Random training set.
    TrainingSet = {'685058', '734045', '826353', '887373', '100408', '110411', '126325', '142828'};
elseif(strcmp(tsid, 'Random3'))
    % Random training set.
    TrainingSet = {'151526', '159239', '169343', '197550', '208226', '217126', '371843', '530635'};
elseif(strcmp(tsid, 'Random4'))
    % Random training set.
    TrainingSet = {'598568', '688569', '856766', '959574', '118730', '127933', '134324', '143325'};
elseif(strcmp(tsid, 'Random5'))
    % Random training set.
    TrainingSet = {'151627', '172332', '190031', '198451', '217429', '284646', '541943', '627549'};
elseif(strcmp(tsid, 'Random6'))
    % Random training set.
    TrainingSet = {'702133', '894673', '978578', '102816', '111312', '118932', '128632', '135932'};
elseif(strcmp(tsid, 'Random7'))
    % Random training set.
    TrainingSet = {'144226', '175439', '191437', '199150', '209935', '221319', '293748', '397760'};
elseif(strcmp(tsid, 'Random8'))
    % Random training set.
    TrainingSet = {'638049', '704238', '751348', '859671', '896879', '984472', '111514', '119833'};
elseif(strcmp(tsid, 'RandomLarge'))
    % Random training set.
    TrainingSet = {'129028', '136833', '148032', '153429', '161731', '176542', '192439', '210617', '224022', '298051', '414229', '547046', '645551', '753251', '861456', '899885',...
                   '103515', '111716', '130013', '137128', '154431', '162329', '177746', '192540', '200614', '239944', '304020', '429040', '559053', '756055', '865363', '901139',...
                   '103818', '112819', '120212', '130316', '138231', '149337', '156233', '193239', '201111', '245333', '307127', '561242', '665254', '715647', '761957', '871964',...
                   '105115', '130922', '138534', '149539', '156637', '182739', '194140', '201414', '212318', '246133', '329440', '485757', '672756', '917255',...
                   '105216', '113619', '123117', '131924', '139637', '150423', '157336', '163432', '195647', '214019', '249947', '497865', '579665', '677968', '729557', '788876', '877168', '932554',...
                   '115320', '124422', '133625', '140420', '150524', '158035', '185139', '196144', '214221', '250427', '355239', '499566', '680957', '732243', '792564', '885975', '937160'};
elseif(strcmp(tsid, 'Monkey'))
    % Diverse training set.
    %TrainingSet = {'Monkey0609'};
    TrainingSet = {'ME3429',  'ME3481',  'ME3487',  'ME3494',  'ME3501',  'ME3513',  'ME3523',  'ME3531',  'ME3537',  'ME3546',  'ME3557',  'ME3583',  'ME3590'};
    ExtraPath = '/';
elseif(strcmp(tsid, 'Lifespan'))
    % Diverse training set.
    TrainingSet = {'LS2001', 'LS2008', 'LS2037', 'LS3017', 'LS3026', 'LS3040', 'LS4025', 'LS4043', 'LS5007', 'LS5040', 'LS5049', 'LS6006', 'LS6038', 'LS2003', 'LS2009', 'LS2043', 'LS3019', 'LS3029', 'LS3046', 'LS4036', 'LS4047', 'LS5038', 'LS5041', 'LS6003', 'LS6009', 'LS6046'};
    ExtraPath = '/Diffusion/Diffusion/';
else
    error('Unknown training set');
end

if(nargin<4)
    ts = length(TrainingSet);
end

if(nargin<5)
    rnds = 4;
end

if(nargin<6)
    ds = 2;
end

if(nargin<7)
    m = ds;
end

DataPath = '/cs/research/vision/hcp/DCA_HCP.2013.3_Proc';

for si = 1:rnds;

for i=1:ts
    load([DataPath '/' TrainingSet{i} ExtraPath sprintf('ipatchlibDS%02i_N%02i', ds, n)]);
    load([DataPath '/' TrainingSet{i} ExtraPath sprintf('opatchlibDS%02i_N01_M%02i', ds, m)]);
    if(n>1)
        % Need to subsample the opatchlib to match the ipatchlib
        load([DataPath '/' TrainingSet{i} ExtraPath sprintf('indices_to_N1_indicesDS%02i_N%02i', ds, n)]);
        opatchlib = opatchlib(indices_to_N1_indices,:);
    end
        
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


%% Compute linear mapping from ipatch to opatch from the libraries.

% Linear transformation.
train_in = [comipatchlib'; 1E-3*ones(1,size(comipatchlib,1))]';
M = train_in\comopatchlib; %pinv(comipatchlib)*comopatchlib;

% Covariance of estimate
S = comopatchlib - (M'*train_in')';
S = S'*S/size(comipatchlib,1);
iS = inv(S);

% Save some space...
comipatchlib = single(comipatchlib);
comopatchlib = single(comopatchlib);

filename = sprintf('Linear%sDS%02i_%ix%i_%ix%i_TS%i_SRi%03i_%04i.mat', tsid, ds, 2*n+1, 2*n+1, m, m, ts, SampleRate, si)
save([DataPath '/TrainingData/' filename], 'M', 'S', 'iS');
filename = sprintf('PatchLibs%sDS%02i_%ix%i_%ix%i_TS%i_SRi%03i_%04i.mat', tsid, ds, 2*n+1, 2*n+1, m, m, ts, SampleRate, si)
save([DataPath '/TrainingData/' filename], 'comipatchlib', 'comopatchlib', 'patchlibindices', '-v7.3');

%% Loop over randomizations
end

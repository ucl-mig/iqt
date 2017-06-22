function open_matlabpool()
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

close_matlabpool;

poolsize = 0;
try
    myPool = parpool('local');
    poolsize = myPool.NumWorkers;
catch ME
end

if poolsize==0
    try
        matlabpool open;
        poolsize = matlabpool('size');
    catch ME
    end
end

fprintf('Matlab pool with %i threads running.\n', poolsize);

function close_matlabpool()
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

isClosed = false;
try
    poolobj = gcp('nocreate');
    delete(poolobj);
    fprintf('Deleted Matlab pool.\n');
    isClosed = true;
catch ME
end

if ~isClosed
    try
        matlabpool close;
        fprintf('Deleted Matlab pool.\n');
        isClosed = true;    
    catch ME
    end
end

if ~isClosed
    fprintf(1, 'No matlabpools to close.\n');
end
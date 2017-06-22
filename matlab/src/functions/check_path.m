function check_path(path)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

if ~strcmp(path(end), filesep)
    error('[IQT] Incorrect path-string, does not end with a slash: %s',...
        path)
end
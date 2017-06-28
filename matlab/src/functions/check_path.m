function check_path(pathstr)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

if ~(strcmp(pathstr(end), '\') || strcmp(pathstr(end), '/'))
    error('[IQT] Incorrect path-string, does not end with a slash: %s',...
        pathstr)
end
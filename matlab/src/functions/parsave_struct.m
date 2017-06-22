function parsave_struct(filename, varStruct, optstr)
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

    if (nargin < 3)
    	optstr = '';
    end

    save(filename, '-struct', 'varStruct', optstr);
    
end

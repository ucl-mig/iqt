function parsave(filename, var, optstr)
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

    varstr = inputname(2);
    eval(sprintf('%s = var;', varstr));
    

    save(filename, varstr, optstr);
    
end

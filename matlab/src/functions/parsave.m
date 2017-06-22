function ag_parsave(filename, var, optstr)

    if (nargin < 3)
    	optstr = '';
    end

    varstr = inputname(2);
    eval(sprintf('%s = var;', varstr));
    

    save(filename, varstr, optstr);
    
end

function ag_parsave_varStruct(filename, varStruct, optstr)

    if (nargin < 3)
    	optstr = '';
    end

    save(filename, '-struct', 'varStruct', optstr);
    
end

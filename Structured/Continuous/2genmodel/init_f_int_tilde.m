function ret=init_f_int_tilde(funcs) %funcs is an ox1 array, each having a function
    o=length(funcs);
    fints=cell(o,1);
    parfor k=1:o
        fints{k}=funcs(k);
    end
    ret=fints;

end


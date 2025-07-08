function ret=init_genprocess(o,N,d,T)

    %rn no noise in genprocess

    RR_tilde=init_RR_tilde(o,N,d,T);
    F_tilde=init_F_tilde(RR_tilde);
    vars=["RR_tilde", "F_tilde"];
    vals=[{RR_tilde},{F_tilde}];



    ret=dictionary(vars,vals);
end

%use as
% genprocess=init_genprocess(o,N,d,T);
% RR_tilde=genprocess("RR_tilde"){1};




function ret=init_genprocess(o,N,d,T)

    %rn no noise in genprocess

    R_tilde=init_R_tilde(o,N,d,T);
    F_tilde=init_F_tilde(R_tilde);
    vars=["R_tilde", "F_tilde"];
    vals=[{R_tilde},{F_tilde}];



    ret=dictionary(vars,vals);
end

%use as
% genprocess=init_genprocess(o,N,d,T);
% R_tilde=genprocess("R_tilde"){1};




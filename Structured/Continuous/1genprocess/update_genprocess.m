function ret=update_genprocess(F_tilde,F,R_tilde,action,dt,t)
    %noise not here for now
    F_tilde=update_F_tilde(F_tilde,F,R_tilde,action,dt,t);
    R_tilde=update_R_tilde(R_tilde,F_tilde,dt,t);
    vars=["R_tilde", "F_tilde"];
    vals=[{R_tilde},{F_tilde}];

    ret=dictionary(vars,vals);
end

%use as
% genprocess=update_genprocess(o,N,d,T);
% R_tilde=genprocess("R_tilde"){1};
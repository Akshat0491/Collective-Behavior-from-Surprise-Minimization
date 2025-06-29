function ret=init_R_tilde(o,N,d,T) %initialise R_tilde
    R_tilde=cell(o,1);
    for p=1:o
        R_tilde{p}=5*ones(N,d,T);
    end
    ret=R_tilde;
end

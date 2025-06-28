function ret=init_R_tilde(o,N,T) %initialise R_tilde
    R_tilde=cell(o,1);
    for p=1:o
        R_tilde{p}=zeros(N,2,T);
    end
    ret=R_tilde;
end

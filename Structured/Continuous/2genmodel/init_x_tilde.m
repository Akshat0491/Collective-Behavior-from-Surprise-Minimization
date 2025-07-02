function ret=init_x_tilde(N,o,L,dL,T)
    x_tilde=cell(N,1);
    parfor i=1:N
        x_tilde{i}=init_R_tilde(o,L,dL,T);
    end
    ret=x_tilde;
end
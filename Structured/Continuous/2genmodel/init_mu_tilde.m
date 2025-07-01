function ret=init_mu_tilde(N,o,L,dL,T)
    mu_tilde=cell(N,1);
    parfor i=1:N
        mu_tilde{i}=init_R_tilde(o,L,dL,T);
    end
    ret=mu_tilde;
end
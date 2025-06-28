function ret=init_Y_tilde(o,N,L,T)
    Y_tilde=cell(o,1);
    for p=1:o
        Y_tilde_p=cell(N,1);
        parfor i=1:N
            Y_tilde_p{i}=zeros(L,T);
        end
        Y_tilde{p}=Y_tilde_p;
    end
    ret=Y_tilde;
end

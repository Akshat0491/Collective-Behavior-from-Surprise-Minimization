function ret=init_Y_ext_tilde(N,o,L,dL,T)
    Y_ext_tilde=cell(N,1);
    parfor i=1:N
        Y_ext_tilde{i}=cell(o,1);
        for k=1:o
            Y_ext_tilde{i}{k}=zeros(L,dL,T)
        end
    end
    ret=Y_ext_tilde;
end
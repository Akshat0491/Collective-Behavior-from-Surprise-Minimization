function ret=init_Y_tilde(N,o,L,dl,T)
    Y_tilde=cell(N,1);
    parfor i=1:N
        Y_tilde{i}=cell(o,1);
        for k=1:o
            Y_tilde{i}{k}=zeros(L,dl,T)
        end
    end
    ret=Y_tilde;
end
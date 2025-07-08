function ret=init_RR_tilde(o,N,d,T) %initialise RR_tilde
    RR_tilde=cell(o,1);
    for p=1:o
        RR_tilde{p}=zeros(N,d,T);
    end
    ret=RR_tilde;
end

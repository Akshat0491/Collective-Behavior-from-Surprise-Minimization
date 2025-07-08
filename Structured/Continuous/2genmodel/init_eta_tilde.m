function ret=init_eta_tilde(etas,L,dL,T)
    %etas: an array containign the eta value for each order
    o=length(etas);
    ret=cell(o,1);
    parfor i=1:o
        ret{i}=etas(i)*ones(L,dL,T)
    end

end

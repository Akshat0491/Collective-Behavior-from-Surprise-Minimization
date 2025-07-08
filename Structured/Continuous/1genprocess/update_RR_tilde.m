function ret=update_RR_tilde(RR_tilde,F_tilde,dt,t)

    o=length(RR_tilde); % order of the gnn
    temp=RR_tilde;
    parfor i=1:o
        RR_tilde{i}(:,:,t)=temp{i}(:,:,t-1) + (dt .* F_tilde{i}(:,:,t));
    end
    ret=RR_tilde;

end
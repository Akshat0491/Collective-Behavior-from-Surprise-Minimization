function ret=update_R_tilde(R_tilde,F_tilde,dt,t)

    o=length(R_tilde); % order of the gnn
    temp=R_tilde;
    parfor i=1:o
        R_tilde{i}(:,:,t)=temp{i}(:,:,t-1) + (dt .* F_tilde{i}(:,:,t));
    end
    ret=R_tilde;

end
function ret=get_curr_mu(mu_tilde,t)
    N=length(mu_tilde);
    o=length(mu_tilde{1});
    curr=cell(N,1);
    parfor i=1:N
        curr{i}=cell(o,1);
        for k=1:o
        curr{i}{k}=mu_tilde{i}{k}(:,:,t);
        end
    end
    ret=curr;
end
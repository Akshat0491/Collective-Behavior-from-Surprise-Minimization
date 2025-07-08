function ret=update_genprocess(F_tilde,F,RR_tilde,action,dt,t,s)
    %noise not here for now
    F_tilde=update_F_tilde(F_tilde,F,RR_tilde,action,dt,t);
    RR_tilde=update_RR_tilde(RR_tilde,F_tilde,dt,t);



    for j=1:size(RR_tilde{1},2)
        R = RR_tilde{1}(:,j,t);

        R_new = ((R < s & R>-s) .* (R)) + ...
                ((R < s & R<=-s) .* (s * ones(size(R)))) +...
                ((R >= s & R>-s) .* (-s * ones(size(R))));
        RR_tilde{1}(:,j,t) = R_new;
    end

    vars=["RR_tilde", "F_tilde"];
    vals=[{RR_tilde},{F_tilde}];

    ret=dictionary(vars,vals);
end

%use as
% genprocess=update_genprocess(o,N,d,T);
% RR_tilde=genprocess("RR_tilde"){1};
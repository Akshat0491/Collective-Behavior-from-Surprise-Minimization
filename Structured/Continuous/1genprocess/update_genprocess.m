function ret=update_genprocess(F_tilde,F,R_tilde,action,dt,t,s)
    %noise not here for now
    F_tilde=update_F_tilde(F_tilde,F,R_tilde,action,dt,t);
    R_tilde=update_R_tilde(R_tilde,F_tilde,dt,t);



    for j=1:size(R_tilde{1},2)
        R = R_tilde{1}(:,j,t);

        R_new = ((R < s & R>-s) .* (R)) + ...
                ((R < s & R<=-s) .* (s * ones(size(R)))) +...
                ((R >= s & R>-s) .* (-s * ones(size(R))));
        R_tilde{1}(:,j,t) = R_new;
    end

    vars=["R_tilde", "F_tilde"];
    vals=[{R_tilde},{F_tilde}];

    ret=dictionary(vars,vals);
end

%use as
% genprocess=update_genprocess(o,N,d,T);
% R_tilde=genprocess("R_tilde"){1};
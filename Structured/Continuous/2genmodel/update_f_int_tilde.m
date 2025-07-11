function ret=update_f_int_tilde(f_int_tilde,mu_tilde_x_i,mu_tilde_v_i)
    %mu_tilde_i       : ox1 cell, each having LxdL, for ith agent
    %f_int_tilde    : ox1 cell each having a function f_int
    % returns       : ox1 cell
    o=length(mu_tilde_x_i);
    ret=cell(o,1);
    parfor k=1:o
        ret{k}=f_int_tilde{k}(mu_tilde_x_i{k},mu_tilde_v_i{k});
    end

end
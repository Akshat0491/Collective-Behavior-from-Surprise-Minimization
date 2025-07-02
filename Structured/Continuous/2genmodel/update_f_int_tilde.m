function ret=f_int_tilde(f_int_tilde,mu_tilde_x,mu_tilde_v)
    %argument       : ox1 cell
    %f_int_tilde    : ox1 cell each having a function f_int
    % returns       : ox1 cell
    o=length(mu_tilde_x);
    ret=cell(o,1);
    parfor k=1:o
        ret{k}=f_int_tilde{k}(mu_tilde_x{k},mu_tilde_v{k})
    end

end
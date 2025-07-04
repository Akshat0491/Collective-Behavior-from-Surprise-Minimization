function ret=init_genmodel(N,o,L,dL,T,gamma_vecdL_x,lambda_vecdL_x,gamma_vecdL_v,lambda_vecdL_v)
    
    % y_int_tilde =init_Y_tilde(N,o,L,dL,T);
    % x_tilde     =init_x_tilde(N,o,L,dL,T);        is it needed?   
    % v_tilde     =init_v_tilde(N,o,L,dL,T);
    

    mu_tilde_x  =init_mu_tilde(N,o,L,dL,T);
    PI_tilde_x  =init_PI_tilde(o,L,gamma_vecdL_x,lambda_vecdL_x);
    mu_tilde_v  =init_mu_tilde(N,o,L,dL,T);
    PI_tilde_v  =init_PI_tilde(o,L,gamma_vecdL_v,lambda_vecdL_v);
    
    
    vfe         =init_vfe(N,T);
   

    vars = ["mu_tilde_x", "PI_tilde_x", "mu_tilde_v", "PI_tilde_v", "vfe"];
    vals = [{mu_tilde_x}, {PI_tilde_x}, {mu_tilde_v}, {PI_tilde_v}, {vfe}];
    ret=dictionary(vars,vals);
end






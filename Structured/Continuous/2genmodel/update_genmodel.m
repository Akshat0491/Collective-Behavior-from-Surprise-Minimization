function ret=update_genmodel()
    
    % y_int_tilde =update_Y_tilde(N,o,L,dL,T);
    % x_tilde     =update_x_tilde(N,o,L,dL,T);        is it needed?   
    % v_tilde     =update_v_tilde(N,o,L,dL,T);
    

    mu_tilde_x  =update_mu_tilde();   %not ready yet, need to implement grad descent
    PI_x        =update_PI_tilde();   %not ready yet, need to implement grad descent
    mu_tilde_v  =update_mu_tilde();   %not ready yet, need to implement grad descent
    PI_v        =update_PI_tilde();   %not ready yet, need to implement grad descent
    
    
    vfe         =update_vfe(vfe,Y_ext_tilde,PI_tilde_y,mu_tilde_x,PI_tilde_x,mu_tilde_v,PI_tilde_v,g_int_tilde,f_int_tilde,eta_tilde,t);
   

    vars = ["mu_tilde_x", "PI_x", "mu_tilde_v", "PI_v", "vfe"];
    vals = [{mu_tilde_x}, {PI_x}, {mu_tilde_v}, {PI_v}, {vfe}];
    ret=dictionary(vars,vals);
end
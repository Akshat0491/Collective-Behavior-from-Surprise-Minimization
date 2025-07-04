% function ret=update_genmodel(vfe,Y_ext_tilde,PI_tilde_y,mu_tilde_x,kx,PI_tilde_x,kpx,mu_tilde_v,kv,PI_tilde_v,kpv,g_int_tilde,jacobian_g_int_tilde_x, jacobian_g_int_tilde_v,f_int_tilde,jacobian_f_int_tilde_x,jacobian_f_int_tilde_v,eta_tilde,t)
    
%     % y_int_tilde =update_Y_tilde();
%     % x_tilde     =update_x_tilde();        is it needed?   
%     % v_tilde     =update_v_tilde();
    

%     mu_tilde_x  =update_mu_tilde_x(mu_tilde_x,kx);   %not ready yet, need to implement grad descent
%     PI_tilde_x  =update_PI_tilde_x(PI_tilde_x,kpx);   %not ready yet, need to implement grad descent
%     mu_tilde_v  =update_mu_tilde_v(mu_tilde_v,kv);   %not ready yet, need to implement grad descent
%     PI_tilde_v  =update_PI_tilde_v(PI_tilde_v,kpv);   %not ready yet, need to implement grad descent
    
    
%     vfe         =update_vfe(vfe,Y_ext_tilde,PI_tilde_y,mu_tilde_x,PI_tilde_x,mu_tilde_v,PI_tilde_v,g_int_tilde,f_int_tilde,eta_tilde,t);
   

%     vars = ["mu_tilde_x", "PI_tilde_x", "mu_tilde_v", "PI_tilde_v", "vfe"];
%     vals = [{mu_tilde_x}, {PI_tilde_x}, {mu_tilde_v}, {PI_v}, {vfe}];
%     ret=dictionary(vars,vals);
% end
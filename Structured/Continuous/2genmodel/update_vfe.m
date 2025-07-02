function ret=update_vfe(vfe,Y_ext_tilde,PI_tilde_y,mu_tilde_x,PI_tilde_x,mu_tilde_v,PI_tilde_v,g_int_tilde,f_int_tilde,eta_tilde,t)
   %g_int_tilde: maps a cell(o,1) to cell(o,1)

    o=length(Y_ext_tilde);
    N=length(mu_tilde_x);
    % mx=mu_tilde_x;
    parfor i=1:N
        
        % e_tilde_y=cell(o,1);
        % e_tilde_x=cell(o,1);
        % e_tilde_v=cell(o,1);
        
        % for k=1:o   
        % e_tilde_y{k}=Y_ext_tilde{i}{k}(:,:,t)-g_int_tilde(mu_tilde_x{i}{k}(:,:,t),mu_tilde_v{i}{k}(:,:,t)); %considering g: (LxdL)x(LxdL) :--> (LxdL);
        % e_tilde_x{k}=D(mu_tilde_x){k}-f_int_tilde(mu_tilde_x{i}{k},mu_tilde_v{i}{k})
        % e_tilde_v{k}=
        % end


        e_tilde_y=cell2mat(Y_ext_tilde{i})-cell2mat(update_g_int_tilde(g_int_tilde,mu_tilde_x{i},mu_tilde_v{i}));
        % e_tilde_y=e_tilde_y(:,:,t);
        
        e_tilde_x=cell2mat(D(mu_tilde_x{i}))-cell2mat(update_f_int_tilde(f_int_tilde,mu_tilde_x{i},mu_tilde_v{i})); 
        e_tilde_v=cell2mat(mu_tilde_v{i})-eta_tilde;      
        % e_tilde_v=e_tilde_v(:,:,t); 



        vfe_it= 0.5*(  a_dot_PI_b(e_tilde_x,PI_tilde_x,e_tilde_x,t)+...
                       a_dot_PI_b(e_tilde_v,PI_tilde_v,e_tilde_v,t)+...
                       a_dot_PI_b(e_tilde_y,PI_tilde_y,e_tilde_y,t));
        vfe(i,t)=sum(vfe_it,2);
    end

    ret=vfe;
end











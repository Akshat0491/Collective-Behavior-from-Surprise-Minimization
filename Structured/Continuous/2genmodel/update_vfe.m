function ret=update_vfe(vfe,Y_ext_tilde,PI_tilde_y,mu_tilde_x,PI_tilde_x,mu_tilde_v,PI_tilde_v,g_int_tilde,f_int_tilde,eta_tilde,t)
   %g_int_tilde: maps a cell(o,1) to cell(o,1)
    dL=1;
    N=length(mu_tilde_x);
    L=size(mu_tilde_x{1}{1},1);
    o=length(mu_tilde_x{1});


    % mx=mu_tilde_x;

    curr_Y_ext_tilde=get_curr_sense(Y_ext_tilde,t);
    curr_mu_tilde_x =get_curr_mu(mu_tilde_x,t-1);
    D_curr_mu_tilde_x =D(get_curr_mu(mu_tilde_x,t-1));
    curr_mu_tilde_v =get_curr_mu(mu_tilde_v,t-1);

    parfor i=1:N
        
        % e_tilde_y=cell(o,1);
        % e_tilde_x=cell(o,1);
        % e_tilde_v=cell(o,1);
        
        % for k=1:o   
        % e_tilde_y{k}=Y_ext_tilde{i}{k}(:,:,t)-g_int_tilde(mu_tilde_x{i}{k}(:,:,t),mu_tilde_v{i}{k}(:,:,t)); %considering g: (LxdL)x(LxdL) :--> (LxdL);
        % e_tilde_x{k}=D(mu_tilde_x){k}-f_int_tilde(mu_tilde_x{i}{k},mu_tilde_v{i}{k})
        % e_tilde_v{k}=
        % end

        gx=cell2mat(update_g_int_tilde(g_int_tilde,curr_mu_tilde_x{i},curr_mu_tilde_v{i}));
        fx=cell2mat(update_f_int_tilde(f_int_tilde,curr_mu_tilde_x{i},curr_mu_tilde_v{i}));

        e_tilde_y=mat2cell(cell2mat(curr_Y_ext_tilde{i})-gx,L*ones(1,o),ones(1,dL));        
        e_tilde_x=mat2cell(cell2mat(D(curr_mu_tilde_x{i}))-fx,L*ones(1,o),ones(1,dL));
        e_tilde_v=cell2mat(curr_mu_tilde_v{i})-eta_tilde; %need to figure this out, also update in update_mu_tilde functions


        vfe_it= 0.5*(  a_dot_PI_b(e_tilde_x,PI_tilde_x,e_tilde_x)+...
                       a_dot_PI_b(e_tilde_v,PI_tilde_v,e_tilde_v)+...
                       a_dot_PI_b(e_tilde_y,PI_tilde_y,e_tilde_y));
        vfe(i,t)=sum(vfe_it,2);
    end

    ret=vfe;
end











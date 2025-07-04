function ret=update_genmodel(vfe,mu_tilde_x,mu_tilde_v,jacobian_g_int_tilde_x,PI_tilde_y,Y_ext_tilde,g_int_tilde,PI_tilde_x,jacobian_f_int_tilde_x,f_int_tilde,eta_tilde,PI_tilde_v,t,dt,kx)
    % flag: string indicating whether it's for x or v
    N=length(mu_tilde_x);
    L=size(mu_tilde_x{1}{1},1);
    dL=size(mu_tilde_x{1}{1},2);
    o=length(mu_tilde_x{1});

    curr_Y_ext_tilde=get_curr_sense(Y_ext_tilde,t);
    curr_mu_tilde_x =get_curr_mu(mu_tilde_x,t-1);
    
    curr_mu_tilde_v =get_curr_mu(mu_tilde_v,t-1);
    % D_curr_mu_tilde_v =D(get_curr_mu(mu_tilde_v,t-1));
    
    parfor i=1:N
    % for i=1:N
        i,t
        gx=cell2mat(update_g_int_tilde(g_int_tilde,curr_mu_tilde_x{i},curr_mu_tilde_v{i}));
        fx=cell2mat(update_f_int_tilde(f_int_tilde,curr_mu_tilde_x{i},curr_mu_tilde_v{i}));

        e_tilde_y=mat2cell(cell2mat(curr_Y_ext_tilde{i})-gx,L*ones(1,o),ones(1,dL));        
        e_tilde_x=mat2cell(cell2mat(D(curr_mu_tilde_x{i}))-fx,L*ones(1,o),ones(1,dL));
        % e_tilde_v=cell2mat(curr_mu_tilde_v{i})-eta_tilde;  
    
        jacobi_gx=make_jacobi_matrix(o,L,jacobian_g_int_tilde_x,gx);
        jacobi_fx=make_jacobi_matrix(o,L,jacobian_f_int_tilde_x,fx);


        grad_term_x_1=jacobi_dot_PI_b(jacobi_gx,PI_tilde_y,e_tilde_y); 
        grad_term_x_2=D_dot_PI_e(PI_tilde_x,e_tilde_x,o); %Keep in mind that this term should be subtracted
        grad_term_x_3=jacobi_dot_PI_b(jacobi_fx,PI_tilde_x,e_tilde_x);

        % grad_term_v_1=jacobi_dot_PI_b(jacobi_gx,PI_tilde_y,e_tilde_y) 
        % grad_term_v_2=D_dot_PI_e(PI_tilde_x,e_tilde_x,o) %Keep in mind that this term should be subtracted
        % grad_term_v_3=jacobi_dot_PI_b(jacobi_fx,PI_tilde_x,e_tilde_x)
        % grad_terms_v=jacobi_dot_PI_b(jacobi_gv,PI_tilde_y,e_tilde_y) - jacobi_dot_PI_b(jacobi_fv,PI_tilde_x,e_tilde_x) + PI_tilde_v{1}*e_tilde_v; %last term is manual currently
        D_curr_mu_tilde_x =D(curr_mu_tilde_x{i});
        for k=1:o
            mu_tilde_x{i}{k}(:,:,t)=curr_mu_tilde_x{i}{k}(:,:) + kx .* (D_curr_mu_tilde_x{k}(:,:) + ( grad_term_x_1{k}(:,:) - grad_term_x_2{k}(:,:) + grad_term_x_3{k}(:,:) ));
            % mu_tilde_v{i}{k}(:,:,t)=curr_mu_tilde_v{i}{k}(:,:); %+ 0 * (D_curr_mu_tilde_v{i}{k}(:,:) + ( grad_term_v_1{k}(:,:) - grad_term_v_2{k}(:,:) + grad_term_v_3{k}(:,:) ));
        end

        vfe_it= 0.5*(  a_dot_PI_b(e_tilde_x,PI_tilde_x,e_tilde_x)+...
                       a_dot_PI_b(e_tilde_y,PI_tilde_y,e_tilde_y)); %+... a_dot_PI_b(e_tilde_v,PI_tilde_v,e_tilde_v,t)
        vfe(i,t)=sum(vfe_it,2);

    end
    vars = ["mu_tilde_x", "PI_tilde_x", "mu_tilde_v", "PI_tilde_v", "vfe"];
    vals = [{mu_tilde_x}, {PI_tilde_x}, {mu_tilde_v}, {PI_tilde_v}, {vfe}];
    ret=dictionary(vars,vals);
end














































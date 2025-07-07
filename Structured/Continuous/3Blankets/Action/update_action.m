function ret=update_action(R_tilde,action,mu_tilde_x,mu_tilde_v,g_int_tilde,Y_ext_tilde,jacobian_G_ext_tilde,PI_tilde_y,ka,dt,t)
    %act_dot will come from VFE/somewhere
    %act_dot : Nxd matrix, where N is the number of agents and d is the dimension of the world
    % action(:,:,t)=action(:,:,t-1)+(dt .* act_dot);
    % action(:,:,t)=action(:,:,t-1);
    % action(:,:,t)= -(ka/dt)*grad_object(F); %something of this form only
    
    % N=size(R_tilde{1},1);
    % action(:,:,t)= zeros(size(action(:,:,1)));


    %playing in 1D from here on

    % d_1_2=R_tilde{1}(1,1,t)-R_tilde{1}(2,1,t)+eps;
    % d_1_2_cap=(d_1_2*1e-31) ./ sqrt(sum(d_1_2 .^2,2));


    % ret=action;
















    N=size(R_tilde{1},1);
    L=size(mu_tilde_x{1}{1},1);
    dL=size(mu_tilde_x{1}{1},2);
    o=length(R_tilde);

    curr_Y_ext_tilde=get_curr_sense(Y_ext_tilde,t-1);
    curr_mu_tilde_x =get_curr_mu(mu_tilde_x,t-1);
    
    curr_mu_tilde_v =get_curr_mu(mu_tilde_v,t-1);

    for i=1:N
    gx=cell2mat(update_g_int_tilde(g_int_tilde,curr_mu_tilde_x{i},curr_mu_tilde_v{i}));
    Y=cell2mat(curr_Y_ext_tilde{i});
    e_tilde_y=mat2cell(Y-gx,L*ones(1,o),ones(1,dL));  

    jacobi_Gext=make_jacobi_matrix_action(i,t,o,L,jacobian_G_ext_tilde,Y,R_tilde);

    acti=jacobi_dot_PI_b(jacobi_Gext,PI_tilde_y,e_tilde_y);
    
        for d=1:size(action,2)
        action(i,d,t)=(i==1).*(ka/dt)*sum(acti{d},1);
        end
    end


    dist_x_21 = R_tilde{1}(2,1,t-1)-R_tilde{1}(1,1,t-1)+eps;

    action(1,1,t)=action(1,1,t)+ (1-R_tilde{2}(1,1,t-1)) + (dist_x_21<8) .* (dist_x_21*-10) ./ norm(dist_x_21)^3;  %agent 1 prefers the velocity of 1 i_cap
    % % % action(1,:,t)= - R_tilde{1}(1,:,t-1);
    % action(2,1,t)=action(2,1,t)+ (-1-R_tilde{2}(2,1,t-1))+ (dist_x_21<8) .* (dist_x_21*10) ./ norm(dist_x_21)^3; %agent 2 prefers the velocity of -1 i_cap


ret=action;




end








function ret=make_jacobi_matrix_action(i,t,o,L,jacobian_func_int_tilde,vals,R_tilde)
    %vals : o*L x1 array. L;L;L...
            vals=reshape(vals,[L,o])';   %makes in oxL
            sz=size(jacobian_func_int_tilde);
        jacobi_mat=cell(sz);    %each cell shal get Lx1 thing
        for m=1:sz(1)
            for n=1:sz(2)
                if isscalar(jacobian_func_int_tilde{m,n}(R_tilde{m}(:,1,t-1)+eps,R_tilde{m}(i,1,t-1),eps+R_tilde{m}(:,2,t-1),R_tilde{m}(i,2,t-1))')
                    jacobi_mat{m,n}=jacobian_func_int_tilde{m,n}(R_tilde{m}(:,1,t-1)+eps,R_tilde{m}(i,1,t-1),eps+R_tilde{m}(:,2,t-1),R_tilde{m}(i,2,t-1))*ones(1,L)';
                else
        jacobi_mat{m,n}=jacobian_func_int_tilde{m,n}(R_tilde{m}(:,1,t-1)+eps,R_tilde{m}(i,1,t-1),eps+R_tilde{m}(:,2,t-1),R_tilde{m}(i,2,t-1));
                end
            end
        end
        ret=jacobi_mat;
    end
    















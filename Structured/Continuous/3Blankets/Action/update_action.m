function ret=update_action(R_tilde,action_obj,t)
    %act_dot will come from VFE/somewhere
    %act_dot : Nxd matrix, where N is the number of agents and d is the dimension of the world
    % action_obj(:,:,t)=action_obj(:,:,t-1)+(dt .* act_dot);
    % action_obj(:,:,t)=action_obj(:,:,t-1);
    % action_obj(:,:,t)= -(ka/dt)*grad_object(F); %something of this form only
    
    action_obj(:,:,t)= 10*zeros(size(action_obj(:,:,1)));
    action_obj(1,:,t)= (t<4).*ones(size(action_obj(1,:,1))) + (t>15) .*(-1)*ones(size(action_obj(1,:,1))) + (t>20) .* ones(size(action_obj(1,:,1))); 
    ret=action_obj;
end















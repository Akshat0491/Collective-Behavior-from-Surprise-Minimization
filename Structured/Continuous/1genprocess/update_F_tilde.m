function ret=update_F_tilde(F_tilde,F,R_tilde,action,dt,t)
    % F : function handle that takes R_tilde and action as inputs and returns the next state of the world.
    % F_tilde : cell array of size o, where o is the order of the gnn.
    % F_tilde{1} : NxdxT matrix, where N is the number of agents and T is the number of time steps.
    % F must be the function for how the last term in gnn approximation, evolves, for example acceleration when x_tilde={x,xdot,xdotdot}
    % R_tilde : cell array of size o, where o is the order of the gnn.
    % R_tilde{1} : NxdxT matrix, where N is the number of agents and T is the number of time steps.
    % R_tilde is taken to be 2D here. Can modify to make 3D or even 1D
    % action=last entry of gnn. Nxd
    % returns a cell array having world's states and the observations
    % R_tilde initialised cell{o,1} : ox1 cell-matrix, each cell containing Nx2xTmatric, where o is the order of gnn, N is the number of agents and T is the number of time steps.
    
    % action : Nx2xT matrix, where N is the number of agents and T is the number of time steps.


    
    o=length(R_tilde); % order of the gnn
    sz=size(R_tilde{1}(:,:,1));
    

    for i=1:o % should be successive
        if i==1                                             
            F_tilde{o-i+1}(:,:,t)= zeros(sz);           %acknowledgement: this line was debuged using chatGPT
        elseif i==2
            F_tilde{o-i+1}(:,:,t)= F(R_tilde,action(:,:,t));
            % F_tilde{o-i+1}(:,:,t)= F(R_tilde,action);
        else
            F_tilde{o-i+1}(:,:,t)= F_tilde{o-i+1}(:,:,t-1) + (dt .* F_tilde{o-i+2}(:,:,t));
        end
    end

    ret=F_tilde;
end



































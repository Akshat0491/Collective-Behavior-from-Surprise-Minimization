function ret=genprocess_gnn(F,R_tilde,action,noise_params,t,dt)  
    % Decide whether to have noise or not on case by case basis
    % F : function handle that takes R_tilde and action as inputs and returns the next state of the world.
    % F must be the function for how the last term in gnn approximation, evolves, for example accleration when x_tilde={x,xdot,xdotdot}
    % R_tilde : cell array of size o, where o is the order of the gnn.
    % R_tilde{1} : Nx2xT matrix, where N is the number of agents and T is the number of time steps.
    % R_tilde is taken to be 2D here. Can modify to make 3D or even 1D
    %action=last entry of gnn. Nx2
    %returns a cell array having world's states and the observations
    % R_tilde initialised cell{o,1} : ox1 cell-matrix, each cell containing Nx2xTmatric, where o is the order of gnn, N is the number of agents and T is the number of time steps.
    % action : Nx2xT matrix, where N is the number of agents and T is the number of time steps.
    % t : current time step.
    % default case, same as book
    %X is the position vectors of every agent.


    %add noises as and when required.
    o=length(R_tilde); %order of the gnn
    sz=size(R_tilde{1});
    F=@(X,action) action;         %should return Nx1 this function governs how the world evolves and the level of last term of R_tilde
                    %this functions maps the world state to the observation


    F_tilde=cell(o,1);
    %initialising F_tilde
    parfor i=1:o
        F_tilde{i}=zeros(sz);
    end

    for i=1:o %should be successive
        if i==1
                F_tilde{o}(:,:,t+1)= F(R_tilde,action);
        else
            F_tilde{o-i+1}=F_tilde{o-i+1} + (dt .* F_tilde{o-i+2});
        end
        R_tilde{o-i+1}(:,:,t+1)=R_tilde{o-i+1}(:,:,t) + (dt .* F_tilde{o}(:,:,t+1));


        

    end


    ret={R_tilde,F_tilde};
end



%define F carefulle, it must contain details of how velocities onwards changes happens, it's last row must come from 
% active inference, and rest from integration@







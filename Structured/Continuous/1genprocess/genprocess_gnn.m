function ret=genprocess_gnn(x_tilde,action,t,dt) 
    %action=last entry of gnn
    %returns a cell array having world's states and the observations
    % X initialised cell{o,1} : ox1 cell-matrix, each cell containing NxTmatric, where o is the order of gnn, N is the number of agents and T is the number of time steps.
    % action : NxT matrix, where N is the number of agents and T is the number of time steps.
    % t : current time step.
    % default case, same as book
    %X is the position vectors of every agent.

    %add noises as and when required.
    o=length(x_tilde); %order of the gnn
    sz=size(x_tilde{1});
    F=@(X,action) action;         %should return Nx1 this function governs how the world evolves and the level of last term of x_tilde
    G_tilde=@(X) X;               %this functions maps the world state to the observation


    F_tilde=cell(o,1);
    %initialising F_tilde
    parfor i=1:o
        F_tilde{i}=zeros(sz);
    end

    for i=1:o %should be successive
        if i==1
                F_tilde{o}(:,t+1)= F(x_tilde,action);
        else
            F_tilde{o-i+1}=F_tilde{o-i+1} + (dt .* F_tilde{o-i+2});
        end
        x_tilde{o-i+1}(:,t+1)=x_tilde{o-i+1}(:,t) + (dt .* F_tilde{o}(:,t+1));


        

    end


    y_tilde=G_tilde(x_tilde);

    ret={x_tilde,y_tilde,F_tilde};
end



%define F carefulle, it must contain details of how velocities onwards changes happens, it's last row must come from 
% active inference, and rest from integration@







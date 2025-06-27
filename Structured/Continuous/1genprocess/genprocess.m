function ret=genprocess(X,action,t,dt) %returns a cell array having world's states and the observations
    % X : NxT matrix, where N is the number of agents and T is the number of time steps.
    % action : NxT matrix, where N is the number of agents and T is the number of time steps.
    % t : current time step.
    % default case, same as book
    %X is the position vectors of every agent.

    %add noises as and when required.
    F=@(X,action) action;  %this function governs how the world evolves
    G=@(X) X;               %this functions maps the world state to the observation

    X(:,t+1)= X(:,t)+ (dt * F(X(:,t),action(:,t)));  %update world dynamics
    % X(:,t+1)=X(:,t)+integral(squeeze(F(X(:,t),action(:,t))),0,dt); Didn't work, but would be cool if it worked. %integral of the function F over the time period dt, uses integrate_T^T+dt=integral_0^dt 
    Y=G(X);

    ret={X,Y};
end





function ret=genmodel(X,action,t) %returns a cell array having world's states and the observations
    % default case, same as book

    %add noises as and when required.
    f=@(X,actiont) action;  %this function governs how the world evolves
    g=@(X) X;               %this functions maps the world state to the observation

    X(t+1)= X(t)+ (dt * f(X(t),action(t)));  %update world dynamics

    Y=G(X);

    ret={X,Y};
end

% f(mue,v), v would be a parameter



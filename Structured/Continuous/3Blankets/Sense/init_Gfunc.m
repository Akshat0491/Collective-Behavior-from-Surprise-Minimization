function ret=init_Gfunc
    % A is a Nx2 matrix, where N is the number of agents
    % G maps the position vector to the observation
    % For example, if A is a position vector, G could be the distance from the origin.
    
    % Example implementation: return the Euclidean norm of each agent's position
    ret = @(A,i,L) sqrt(sum((A(1:L)-A(i)).^2, 2));  % returns a column vector of size Lx1, 
end




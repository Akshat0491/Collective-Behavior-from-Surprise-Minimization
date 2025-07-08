function ret=init_F_tilde(RR_tilde)
    % F : function handle that takes RR_tilde and action as inputs and returns the next state of the world.
    % F must be the function for how the last term in gnn approximation, evolves, for example acceleration when x_tilde={x,xdot,xdotdot}
    % RR_tilde : cell array of size o, where o is the order of the gnn.
    % RR_tilde{1} : Nx2xT matrix, where N is the number of agents and T is the number of time steps.
    % RR_tilde is taken to be 2D here. Can modify to make 3D or even 1D
    % action=last entry of gnn. Nx2
    % returns a cell array having world's states and the observations
    % RR_tilde initialised cell{o,1} : ox1 cell-matrix, each cell containing Nx2xTmatric, where o is the order of gnn, N is the number of agents and T is the number of time steps.
    % action : Nx2xT matrix, where N is the number of agents and T is the number of time steps.

    o=length(RR_tilde); % order of the gnn
    sz=size(RR_tilde{1});

    F_tilde=cell(o,1);
    % initialising F_tilde
    parfor i=1:o
        F_tilde{i}=zeros(sz);
    end
    ret=F_tilde;

end
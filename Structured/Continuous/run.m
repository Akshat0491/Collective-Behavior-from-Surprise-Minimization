o=3;
N=10;
T=100;
dimworld=2; %aka d % dimension of the variable that charecterises the world state.
dt=0.1; % time step
R_tilde=init_R_tilde(o,N,T);

F=@(R_tilde,action) action; % function that governs the last term in the gnn approximation. 
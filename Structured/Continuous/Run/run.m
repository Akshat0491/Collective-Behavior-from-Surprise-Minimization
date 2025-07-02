o=3;
N=10;
T=100;
dimworld=2; %aka d % dimension of the variable that charecterises the world state.
dt=0.1; % time step
R_tilde=init_R_tilde(o,N,T);
L=4; 
dL = 1; %Do not change this for now, might give weird stuff for dL>1, for such a case, use multiple genmodels insetad.


F=@(R_tilde,action) action; % function that governs the last term in the gnn approximation. 

% Example implementation: return the Euclidean norm of each agent's position
G=@(R_tilde,i,L,dL) sqrt((R_tilde(1:L,1:dL)-A(i,1:dL)*0).^2);  % maps Nxd-->LxdL 
%i: ith agent,L: number of informations stored, like sensory sectors of heins et al 2024, dL: dimenstion of the world as agent sees it. For example, agent might be seeing only distances, making the world 1D for it, though the world maybe dD.


g=@(mu_tilde_x,mu_tilde_v) mu_tilde_x;
f=@(mu_tilde_x,v_tilde) 1;
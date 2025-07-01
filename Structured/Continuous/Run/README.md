This folder contains single matlab scripts for differrent differrent contexts. for every differrent context, make a copy of the template


Initialisations/Parameters









o   : Order Upto which the Generalised coordinates are truncated to. Generally o=3
N   : Number of Agents
d   : dimension of the variable that charecterises the world state. Like d=2 for a 2D world (x,y). d=3 for 3D world(x,y,z). 
L   : Number of info pieaces/agents as per individual agent
dl  : dimension of the world as per the agent.

T   : Total Number of Time Steps
dt  : time step

F   : (Nxd)xaction --> (Nxd) Function handle that governs how action translates to 
      last term of the gnn truncation. For exampled if o=3, action=acceleration, and F would be the net force. DR=F(R,action) + noise. Usually: in a world governed only by actions, F=@(R_tilde,action) action; % function that governs the last term in the gnn approximation.....
      arguments. F(R_tilde,action)


G   : (N,d)-->(L,dL). A function handle that converts world information to agents setup. 
       Y=G(R_tilde)+noise. is w.r.t ith agent. Arguments. G(R_tilde,i,L,dl). 

       How ith agent sees R_tilde and fits into Lxdl


Note: No noise is added right now. Until this message is deleted






















Description of Objects involved:

R_tilde : Object containing the worlds information. ox1 cell, each cell having NxdxT matrix describing the world.

mu_tilde : Object containing beliefs about the world.
            Nx1 cell, each cell having ox1 cell, each having LxdLxT matric describing the beliefs about the world.

action: NxdxT

perception: LxdLxT




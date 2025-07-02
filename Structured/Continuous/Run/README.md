This folder contains single matlab scripts for differrent differrent contexts. for every differrent context, make a copy of the template


Initialisations/Parameters









o   : Order Upto which the Generalised coordinates are truncated to. Generally o=3
N   : Number of Agents
d   : dimension of the variable that charecterises the world state. Like d=2 for a 2D world (x,y). d=3 for 3D world(x,y,z). 
L   : Number of info pieaces/agents as per individual agent
dl  : dimension of the world as per the agent. rn defaulting it to 1 everywhere, non-1 implementation is currently uncertain.

T   : Total Number of Time Steps
dt  : time step

gamme_vecdL_x  : gamma  for x for each dimensions dL. [1xdL]
lambda_vecdL_x : lambda for x for each dimensions dL. [1xdL]
gamme_vecdL_v  : gamma  for v for each dimensions dL. [1xdL]
lambda_vecdL_v : lambda for v for each dimensions dL. [1xdL]


F   : ((Nxd),action (1xd)) --> (Nxd) Function handle that governs how action translates to 
      last term of the gnn truncation. For exampled if o=3, action=acceleration, and F would be the net force. DR=F(R,action) + noise. Usually: in a world governed only by actions, F=@(R_tilde,action) action; % function that governs the last term in the gnn approximation.....
      arguments. F(R_tilde,action)


G   : (N,d)-->(L,dL). A function handle that converts world information to agents setup. 
       Y=G(R_tilde)+noise. is w.r.t ith agent. Arguments. G(R_tilde,i,L,dl). 

       How ith agent sees R_tilde and fits into Lxdl


Note: No noise is added right now. Until this message is deleted


f=[] %ox1 array containing functions handles. each function handles maps ((LxdL)x(LxdL))-->LxdL; %doubt here, would this be ((LxdL)x(LxdL))-->Nxd
g=[] %ox1 array containing functions handles. each function handles maps ((LxdL)x(LxdL))-->LxdL;



f_int_tilde       :((ox1 cell)x(ox1 cell))-->(ox1 cell)
g_ini_tilde       : (L,dL)-->(L,dL). A function handle that converts beliefs into predictions. y


be careful that you pass mu_tilde for i'th agent while updating






















Description of Objects involved:

R_tilde : Object containing the worlds information. ox1 cell, each cell having NxdxT matrix describing the world.

mu_tilde : Object containing beliefs about the world.
            Nx1 cell, each cell having ox1 cell, each having LxdLxT matric describing the beliefs about the world.

action: NxdxT

perception: LxdLxT




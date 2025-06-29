% function ret=test_genprocess
%testing a spring maxx system

%parameters
o=3; %truncate upto acceleration
N=2; %Number of agents
d=1; %One dimensional world 
T=200;
dt=0.1;


% F=@(A,B) B;
F = @(R_tilde, action) action;  % F = -x(t), i.e. Hookean force


%initialising
time=0:dt:T;
genprocess=init_genprocess(o,N,d,T); %returns a dictionary containing R_tilde and F_tilde like "R_tilde"-->{R_tilde}
R_tilde=get_from_dict(genprocess,"R_tilde"); %get the R_tilde from the dictionary

R_tilde{1}(:,:,:)=5*ones(N,d,T); 
R_tilde{2}(:,:,:)=zeros(N,d,T);    %manually change initial configurations  ||
R_tilde{3}(:,:,:)=zeros(N,d,T);                                      %      ||


F_tilde=get_from_dict(genprocess,"F_tilde"); %get the F_tilde from the dictionary
action=init_action(N,d,T);

%simulation loop
for t=2:T
    R_tilde=get_from_dict(genprocess,"R_tilde");
    action(:,:,t)= - R_tilde{1}(:,:,t-1) - 0.5*R_tilde{2}(:,:,t-1) ;
    genprocess=update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"R_tilde"),action,dt,t);
end

plot(squeeze(R_tilde{1}(1,1,1:T-1)),squeeze(R_tilde{2}(1,1,1:T-1)));
title("Phase Space Rep")
xlabel("X")
ylabel("X_{dot}")
ret=R_tilde;



% end




% update_genprocess(F_tilde,F,R_tilde,action,dt,t)
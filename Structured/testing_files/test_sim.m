% function ret=test_genprocess
%testing a spring maxx system
close all;
clear all;

%---------------------------------------------------------------------------------------------------------------------------------------
%% parameters
% ---------------------------------------------------------------------------------------------------------------------------------------
dt=0.1;
o=3; %truncate upto acceleration
N=4; %Number of agents
d=2; %One dimensional world 
L=N; %number of distinct infos agent holds
dL=1; %dimensions of the dintinct info, never change is for now, doesn't work for dl>1
T=200;

%learning parameters
kx  =0.1;
kpx =0.01;
kv  =0.01; %learnign rates. 
kpv =0.01;
ka  =0.01;

gamma_vecdL_y=ones(1,dL);
gamma_vecdL_x=ones(1,dL); %can have it differrent
gamma_vecdL_v=ones(1,dL);

lambda_vecdL_y=ones(1,dL);
lambda_vecdL_x=ones(1,dL);
lambda_vecdL_v=ones(1,dL);
noise_params=1;



% F=@(A,B) B;
F = @(R_tilde, action) action;  % F = -x(t), i.e. Hookean force
G = @(A,i) sqrt(sum((A - A(i,:)).^2, 2)); %returns euclidean distance 


%will need to change it as per the value of o
% f_int_tilde= @(mu_tilde_x,mu_tilde_v) { mu_tilde_v{1}-mu_tilde_x{1};
%                                         -mu_tilde_x{2};
%                                         -mu_tilde_x{3}};


% this one was kinda toogood to be true


% g_int_tilde= @(mu_tilde_x,mu_tilde_v) {mu_tilde_x{1};     %ox1 cell. ((LxdL)x(LxdL))xo-->(LxdL)xo;
%                              mu_tilde_x{2};
%                              mu_tilde_x{3};
% };



f_int_tilde = {@(x,v)                  x+v;...            %for now, no coupling of x and v is allowed. v should not appear in delf/delx. for example, x*v^2 not lalowed, atleast for current implementation
               @(x_dot,v_dot)          x_dot+v_dot;...                 %this kind of structure is important, atleasft for the way gradient descent is implemented
               @(x_dot_dot,v_dot_dot)  x_dot_dot+v_dot_dot};            % for now, working only with f^[i] is a function only of ith derivatives of x and v.

g_int_tilde = {@(x,v)                  x;
               @(x_dot,v_dot)          x_dot;                 %this kind of structure is important, atleasft for the way gradient descent is implemented
               @(x_dot_dot,v_dot_dot)  x_dot_dot};





%---------------------------------------------------------------------------------------------------------------------------------------
%% Initialising required objects
%---------------------------------------------------------------------------------------------------------------------------------------


time=0:dt:T;


jacobian_f_int_tilde_x=jacob(f_int_tilde,["x","x_dot","x_dot_dot"]);        %do check what this is actually a function of, can make it a .m file later for further clarity, rn all functions would just be functions of the vectore elements w.r.t whose gradient has been sought. 
jacobian_f_int_tilde_v=jacob(f_int_tilde,["v","v_dot","v_dot_dot"]);        % returns an oxo cell, each having the i-jth elements of jacobian as function handle.

jacobian_g_int_tilde_x=jacob(g_int_tilde,["x","x_dot","x_dot_dot"]);
jacobian_g_int_tilde_v=jacob(g_int_tilde,["v","v_dot","v_dot_dot"]);




%each of following objects returns a dictionary containing object by their literal names like "R_tilde"-->{R_tilde}
sense       =init_sense(N,o,L,dL,T,gamma_vecdL_y,lambda_vecdL_y);  
genmodel    =init_genmodel(N,o,L,dL,T,gamma_vecdL_x,lambda_vecdL_x,gamma_vecdL_v,lambda_vecdL_v);
action      =init_action(N,d,T);
genprocess  =init_genprocess(o,N,d,T); 


%% Change initial configuration. Alternatively change the function init_R_tilde

R_tilde=get_from_dict(genprocess,"R_tilde"); %get the R_tilde from the dictionary
S=50;

for i=1:N
    R_tilde{1}(i,:,:)=(2*mod(i,2)-1+i) * ones(1,d,T); 
    R_tilde{2}(i,:,:)=zeros(1,d,T);    %manually change initial configurations  ||
    R_tilde{3}(i,:,:)=zeros(1,d,T);    %      ||
end
genprocess{"R_tilde"}=R_tilde;

F_tilde=get_from_dict(genprocess,"F_tilde"); %get the F_tilde from the dictionary

%---------------------------------------------------------------------------------------------------------------------------------------
%% simulation loop
%---------------------------------------------------------------------------------------------------------------------------------------


for t=2:T

    
    %each of following objects returns a dictionary containing object by their literal names like "R_tilde"-->{R_tilde}
    sense       =update_sense(sense{"Y_ext_tilde"},G,...
                              genprocess{"R_tilde"},noise_params,sense{"PI_tilde_y"},...
                              t);  
    
    genmodel    =update_genmodel(genmodel{"vfe"},genmodel{"mu_tilde_x"},genmodel{"mu_tilde_v"},jacobian_g_int_tilde_x,sense{"PI_tilde_y"},sense{"Y_ext_tilde"},g_int_tilde,genmodel{"PI_tilde_x"},jacobian_f_int_tilde_x,f_int_tilde,1,genmodel{"PI_tilde_v"},t,dt,kx);
    

    
   
    action      =update_action(genprocess{"R_tilde"},action,t);
    
    genprocess  =update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"R_tilde"),action,dt,t); 
















    % R_tilde=get_from_dict(genprocess,"R_tilde");
    % action(:,:,t)= - R_tilde{1}(:,:,t-1) ; %To be replaced with genmodel("action")/updateaction
    % genprocess=update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"R_tilde"),action,dt,t);
    % Y_ext_tilde
end


%---------------------------------------------------------------------------------------------------------------------------------------
%% Plotting and analysing
%---------------------------------------------------------------------------------------------------------------------------------------
R_tilde=genprocess{"R_tilde"};
vfe=genmodel{"vfe"};
mu=genmodel{"mu_tilde_x"};
Y=sense{"Y_ext_tilde"};
figure(1);
scatter(1:T,squeeze(mu{2}{1}(1,1,:)))
hold on;
scatter(1:T,squeeze(Y{2}{1}(1,1,:)))
title("mu vs sense")
xlabel("X")
ylabel("X_{dot}")


figure(2);
scatter(1:T,vfe(1,:))
title("Vfe of 1st agent")
xlabel("t")
ylabel("y")



% S=10;



output_dir ="D:\Projects\Summer 2025\Surprise Minimisation\outputs"; % replace with your desired path

now_dt = datetime('now');
now_dt.TimeZone = 'local';
now=string(now_dt, 'd-MMM-y HH:mm:ss');
now=split(now,":");
now=join(now,"_");
filename = fullfile(output_dir, sprintf('SHM_%s', ...
now));

outputVideo = VideoWriter(filename,'MPEG-4');
outputVideo.FrameRate = 15; % Adjust as needed
% outputVideo.FileFormat='mp4';
open(outputVideo);

fig = figure('Color','w'); % black background for visual clarity
S=max(R_tilde{1});
for t = 1:T-10
clf

scatter(reshape(R_tilde{1}(:,1,t),[N,1]), reshape(R_tilde{1}(:,2,t),[N,1]), 10, 'w', 'filled'); % white dots
ylabel("Y");
xlim([-S S]);
ylim([-S S]);
% quiver(pos(:,1,t), pos(:,2,t), 2*vel(:,1,t)./(sqrt(sum(vel(:,:,t).^2,2))+eps), 2*vel(:,2,t)./(sqrt(sum(vel(:,:,t).^2,2))+eps),"off","Marker",".","ShowArrowHead","on")
% xlabel(sprintf('X — sensed: %.2f, belief: %.2f', ...
%     mean(sense_dists(:,:,t), 'all'), ...
%     mean(mue_dists(:,:,t), 'all')));



% if have_noise
% title(sprintf('Time %.2f sec,N:%.f,\n dt:%.2f, k:%.2f, R:%.2f, vision:full,\n alpha:%f, eta:%.2f, Lambda_w:%.2f,Lambda_z:%.2f, Gamma_w:%.2f, Gamma_z: %.2f', (t-1)*dt,Na,dt,k,R,alpha1,eta,Lambda_w,Lambda_z,Gamma_w,Gamma_z), 'Color', 'w');
% else
% title(sprintf('Without noise Time %.2f sec,N:%.f,\n dt:%.2f, k:%.2f, R:%.2f, vision:full,\n alpha:%f, eta:%.2f, Lambda_w:%.2f,Lambda_z:%.2f, Gamma_w:%.2f, Gamma_z: %.2f', (t-1)*dt,Na,dt,k,R,alpha1,eta,Lambda_w,Lambda_z,Gamma_w,Gamma_z), 'Color', 'w');
% end
set(gca, 'Color', 'k'); % black background for axes
title(sprintf('Time step: %d', t))
drawnow;

% Capture the plot as a frame
frame = getframe(fig);
writeVideo(outputVideo, frame);
end    
    close(outputVideo);






ret=R_tilde;



% end




% update_genprocess(F_tilde,F,R_tilde,action,dt,t)



















%---------------------------------------------------------------------------------------------------------------------------------------
%% Helper Functions
%---------------------------------------------------------------------------------------------------------------------------------------


















































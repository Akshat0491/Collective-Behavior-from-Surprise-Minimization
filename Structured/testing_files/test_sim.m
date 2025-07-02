% function ret=test_genprocess
%testing a spring maxx system
close all;
clear all;
%parameters
o=3; %truncate upto acceleration
N=5; %Number of agents
d=2; %One dimensional world 
L=3; %number of distinct infos agent holds
dL=1; %dimensions of the dintinct info, never change is for now, doesn't work for dl>1
T=100;
gamma_vecdL_y=ones(1,dL);
gamma_vecdL_x=ones(1,dL); %can have it differrent
gamma_vecdL_v=ones(1,dL);

lambda_vecdL_y=ones(1,dL);
lambda_vecdL_x=ones(1,dL);
lambda_vecdL_v=ones(1,dL);
dt=0.1;


% F=@(A,B) B;
F = @(R_tilde, action) action;  % F = -x(t), i.e. Hookean force
% G = something
f_1
f_2
f_3
f=[]; %ox1 array containing functions handles. each function handles maps ((LxdL)x(LxdL))-->LxdL; %doubt here, would this be ((LxdL)x(LxdL))-->Nxd

g_1
g_2
g_3
g=[]; %ox1 array containing functions handles. each function handles maps ((LxdL)x(LxdL))-->LxdL;

f_int_tilde=init_f_int_tilde(f);
g_int_tilde=init_g_int_tilde(g);





%initialising
time=0:dt:T;


%each of following objects returns a dictionary containing object by their literal names like "R_tilde"-->{R_tilde}
sense       =init_sense(N,o,L,dL,T,gamma_vecdL_y,lambda_vecdL_y);  
genmodel    =init_genmodel(N,o,L,dL,T,gamma_vecdL_x,lambda_vecdL_x,gamma_vecdL_v,lambda_vecdL_v);
action      =init_action(N,d,T);
genprocess  =init_genprocess(o,N,d,T); 


%Change initial configuration. Alternatively change the function init_R_tilde

R_tilde=get_from_dict(genprocess,"R_tilde"); %get the R_tilde from the dictionary
S=10;
R_tilde{1}(:,:,:)=5*ones(N,d,T); 
R_tilde{2}(:,:,:)=zeros(N,d,T);    %manually change initial configurations  ||
R_tilde{3}(:,:,:)=zeros(N,d,T);                                      %      ||


F_tilde=get_from_dict(genprocess,"F_tilde"); %get the F_tilde from the dictionary

%simulation loop
for t=2:T

    
    %each of following objects returns a dictionary containing object by their literal names like "R_tilde"-->{R_tilde}
    sense       =update_sense(get_from_dict(sense,"Y_ext_tilde"),G,get_from_dict(genprocess,"R_tilde"),noise_params,t);  
    genmodel    =update_genmodel();
    action      =update_action();
    genprocess  =update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"R_tilde"),action,dt,t); 
















    R_tilde=get_from_dict(genprocess,"R_tilde");
    action(:,:,t)= - R_tilde{1}(:,:,t-1) ; %To be replaced with genmodel("action")/updateaction
    genprocess=update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"R_tilde"),action,dt,t);
    Y_ext_tilde
end

figure(1);
plot(squeeze(R_tilde{1}(1,1,1:T-1)),squeeze(R_tilde{2}(1,1,1:T-1)));
title("Phase Space Rep")
xlabel("X")
ylabel("X_{dot}")


figure(2);
plot(time(1:T-1),squeeze(R_tilde{1}(1,2,1:T-1)));
title("Time evolution of R_tilde{1} for agent 1")
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
ylabel("Y");
xlim([-S S]);
ylim([-S S]);
for t = 1:T-10
clf
scatter(squeeze(R_tilde{1}(1,1,t)), squeeze(R_tilde{1}(1,2,t)), 10, 'w', 'filled'); % white dots
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
drawnow;

% Capture the plot as a frame
frame = getframe(fig);
writeVideo(outputVideo, frame);
end    
    close(outputVideo);






ret=R_tilde;



% end




% update_genprocess(F_tilde,F,R_tilde,action,dt,t)
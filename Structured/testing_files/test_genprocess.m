% function ret=test_genprocess
%testing a spring maxx system

%parameters
o=3; %truncate upto acceleration
N=2; %Number of agents
d=2; %One dimensional world 
T=200;
dt=0.1;


% F=@(A,B) B;
F = @(R_tilde, action) action;  % F = -x(t), i.e. Hookean force
% G = something


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
    % action(:,:,t)= - R_tilde{1}(:,:,t-1) - 0.5*R_tilde{2}(:,:,t-1) ; To be replaced with update_action
    genprocess=update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"R_tilde"),action,dt,t);
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

fig = figure('Color','k'); % black background for visual clarity

for t = 1:T-10
clf
scatter(squeeze(R_tilde{1}(1,1,t)), squeeze(R_tilde{1}(1,2,t)), 10, 'w', 'filled'); % white dots
% quiver(pos(:,1,t), pos(:,2,t), 2*vel(:,1,t)./(sqrt(sum(vel(:,:,t).^2,2))+eps), 2*vel(:,2,t)./(sqrt(sum(vel(:,:,t).^2,2))+eps),"off","Marker",".","ShowArrowHead","on")
% xlabel(sprintf('X — sensed: %.2f, belief: %.2f', ...
%     mean(sense_dists(:,:,t), 'all'), ...
%     mean(mue_dists(:,:,t), 'all')));

ylabel("Y");
xlim([-S S]);
ylim([-S S]);

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
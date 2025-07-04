% function ret=test_genprocess
%testing a spring maxx system
% close all;
clear all;

%---------------------------------------------------------------------------------------------------------------------------------------
%% parameters
% ---------------------------------------------------------------------------------------------------------------------------------------
dt=0.1;
o=3; %truncate upto acceleration
n=2;
N=n^2; %Number of agents
d=2; %One dimensional world 
L=N; %number of distinct infos agent holds
dL=1; %dimensions of the dintinct info, never change is for now, doesn't work for dl>1
T=100;
wanna_save=0;

%learning parameters
n_kx =5;
kx  =0.15;
kpx =0.01;
kv  =0.01; %learnign rates. 
kpv =0.01;
ka  =0.01;

gamma_vecdL_y=ones(1,dL);
gamma_vecdL_x=1*ones(1,dL); %can have it differrent
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



% for i=1:Na
    
% end


di=1/(n-1);
for i=1:N
    % R_tilde{1}(i,:,:)=(2*mod(i,2)-1+i) * ones(1,d,T); 
    R_tilde{1}(i,:,1)=[-1/2 +  (mod(i-1,n))*di,-1/2 + (floor((i-1)/n))*di]; 
    R_tilde{2}(i,:,:)=zeros(1,d,T);    %manually change initial configurations  ||
    R_tilde{3}(i,:,:)=zeros(1,d,T);    %      ||
end
genprocess{"R_tilde"}=R_tilde;

F_tilde=get_from_dict(genprocess,"F_tilde"); %get the F_tilde from the dictionary

%---------------------------------------------------------------------------------------------------------------------------------------
%% simulation loop
%---------------------------------------------------------------------------------------------------------------------------------------


for t=2:T
    t
    
    %each of following objects returns a dictionary containing object by their literal names like "R_tilde"-->{R_tilde}
    action      =update_action(genprocess{"R_tilde"},action,t);

    genprocess  =update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"R_tilde"),...
                                                action,dt,t);
    
    sense       =update_sense(sense{"Y_ext_tilde"},G,...
                              genprocess{"R_tilde"},noise_params,sense{"PI_tilde_y"},...
                              t);  
    
    genmodel    =update_genmodel(genmodel{"vfe"},genmodel{"mu_tilde_x"},genmodel{"mu_tilde_v"},...
                                jacobian_g_int_tilde_x,sense{"PI_tilde_y"},sense{"Y_ext_tilde"},...
                                g_int_tilde,genmodel{"PI_tilde_x"},jacobian_f_int_tilde_x,f_int_tilde,...
                                1,genmodel{"PI_tilde_v"},t,dt,kx,n_kx);
    



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

fig = figure('Color', 'w');
% Split the long title into multiple lines for better visibility
mainTitle = sprintf('dt=%.2f, o=%d, n=%d, N=%d, d=%d, L=%d, dL=%d, T=%d', dt, o, n, N, d, L, dL, T);
paramTitle = sprintf('n_kx=%d, kx=%.2f, kpx=%.2f, kv=%.2f, kpv=%.2f, ka=%.2f', n_kx, kx, kpx, kv, kpv, ka);
tl = tiledlayout(4,1, 'Padding', 'compact', 'TileSpacing', 'compact');
sg = sgtitle(tl, {mainTitle; paramTitle}, 'FontWeight', 'bold', 'FontSize', 12);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
% tl.TitleTop = 10; % Increase space above the title for visibility

% 1st subplot: Prediction error
nexttile;
plot(1:T, squeeze(abs(mu{1}{1}(2,1,:) - Y{1}{1}(2,1,:))), 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]); hold on;
plot(1:T, squeeze(abs(mu{2}{1}(1,1,:) - Y{2}{1}(1,1,:))), 'LineWidth', 1.5, 'Color', [0.12 0.45 0.70]);
hold off;
title(sprintf('|Prediction error Vs Time.| \n Agent i sees the agent j'), 'Color', [0.85 0.33 0.1]);
xlabel("t");
ylabel('|g(\mu) - y|');
legend('Agent 1', 'Agent 2');
grid on;

% 2nd subplot: VFE
nexttile;
plot(1:T, vfe(1,:), 'LineWidth', 1.5, 'Color', [0.13 0.55 0.13]); hold on;
plot(1:T, vfe(2,:), 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]);
hold off;
title("VFE of agent 1 (green) and agent 2 (red)", 'Color', [0.13 0.55 0.13]);
xlabel("t");
ylabel("VFE");
legend('Agent 1', 'Agent 2');
grid on;

% 3rd subplot: Action profile
nexttile;
scatter(1:T, squeeze(action(1,1,:)), 20, [0.12 0.45 0.70], 'filled', 'MarkerEdgeColor', [0.12 0.45 0.70], 'MarkerFaceAlpha', 0.7); hold on;
scatter(1:T, squeeze(action(2,1,:)), 20, [0.85 0.33 0.1], 'filled', 'MarkerEdgeColor', [0.85 0.33 0.1], 'MarkerFaceAlpha', 0.7);
hold off;
title("Action profile of agent 1 (blue) and 2 (red), x coordinate", 'Color', [0.12 0.45 0.70]);
xlabel("time");
ylabel("Action");
legend('Agent 1', 'Agent 2');
grid on;

% 4th subplot: X coordinates of agent 1 and 2
nexttile;
plot(1:T, squeeze(R_tilde{1}(1,1,:)), 'LineWidth', 1.5, 'Color', [0.2 0.2 0.8]); hold on;
plot(1:T, squeeze(R_tilde{1}(2,1,:)), 'LineWidth', 1.5, 'Color', [0.8 0.2 0.2]);
hold off;
title('X coordinate of agent 1 (blue) and 2 (red)', 'Color', [0.2 0.2 0.8]);
xlabel('time');
ylabel('X position');
legend('Agent 1', 'Agent 2');
grid on;


% S=10;


if wanna_save
output_dir ="D:\Projects\Summer 2025\Surprise Minimisation\outputs\structured_testing"; % replace with your desired path

now_dt = datetime('now');
now_dt.TimeZone = 'local';
now=string(now_dt, 'd-MMM-y HH:mm:ss');
now=split(now,":");
now=join(now,"_");
filename = fullfile(output_dir, sprintf('Test_%s', ...
now));


outputVideo = VideoWriter(filename,'MPEG-4');
outputVideo.FrameRate = 15; % Adjust as needed
% outputVideo.FileFormat='mp4';
open(outputVideo);

fig = figure('Color','w'); % black background for visual clarity
% S=max(squeeze(R_tilde{1}(1,1,:)));
S=1;
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
% writeVideo(outputVideo, frame);
end    
    close(outputVideo);
end






ret=R_tilde;



% end




% update_genprocess(F_tilde,F,R_tilde,action,dt,t)



















%---------------------------------------------------------------------------------------------------------------------------------------
%% Helper Functions
%---------------------------------------------------------------------------------------------------------------------------------------



































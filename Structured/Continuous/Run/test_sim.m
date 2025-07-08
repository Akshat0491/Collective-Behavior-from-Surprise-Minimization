% function ret=test_genprocess
%testing a spring maxx system
clc;
close all;
clear all;

%---------------------------------------------------------------------------------------------------------------------------------------
%% parameters
% ---------------------------------------------------------------------------------------------------------------------------------------
dt=0.1;
o=3; %truncate upto acceleration, by default let it be 3, there might be some problems with non 3 values.
n=2;
N=2; %Number of agents
d=2; %One dimensional world 
L=N; %number of distinct infos agent holds
dL=1; %dimensions of the dintinct info, never change is for now, doesn't work for dl>1
T=200;
% Rescale agent diameter to fit nicely in a 2x2 square
S = 1; % Side length of the square (from -S to S)
dia = 1 * (2*S) * 72; % 10% of the square's width, scaled for scatter marker size
wanna_save=1;
S=10;

%learning parameters
n_kx =5*0;
kx0  =0.05*0;
kpx =0.01*0;
kv0  =0.2*1; %learnign rates. 
kpv =0.01*0;
ka0  =0.1*1;

gamma_vecdL_y=ones(1,dL);
gamma_vecdL_x=1*ones(1,dL); %can have it differrent
gamma_vecdL_v=ones(1,dL);

lambda_vecdL_y=ones(1,dL);
lambda_vecdL_x=ones(1,dL);
lambda_vecdL_v=ones(1,dL);
noise_params=1;
etas=[0,10,0]; %ox1 vector



% F=@(A,B) B;
F = @(RR_tilde, action) action;  % F = -x(t), i.e. Hookean force
G = @(A,i) sqrt(sum((A - A(i,:)).^2, 2)); %returns euclidean distance between each agent and the ith agent.





G_ext_tilde={@(x,xi,y,yi,x_dot,xi_dot,y_dot,yi_dot) sqrt(eps + (x-xi).^2 + (y-yi).^2); 
             @(x,xi,y,yi,x_dot,xi_dot,y_dot,yi_dot) ((x-xi)*(x_dot-xi_dot)+(y-yi)*(y_dot-yi_dot)) ./ sqrt((x_dot-xi_dot).^2 + (y_dot-yi_dot).^2)
             @(x,xi,y,yi,x_dot,xi_dot,y_dot,yi_dot) 0};


%will need to change it as per the value of o
% f_int_tilde= @(mu_tilde_x,mu_tilde_v) { mu_tilde_v{1}-mu_tilde_x{1};
%                                         -mu_tilde_x{2};
%                                         -mu_tilde_x{3}};


% this one was kinda toogood to be true


% g_int_tilde= @(mu_tilde_x,mu_tilde_v) {mu_tilde_x{1};     %ox1 cell. ((LxdL)x(LxdL))xo-->(LxdL)xo;
%                              mu_tilde_x{2};
%                              mu_tilde_x{3};
% };



f_int_tilde = {@(x,v)                  x-v;...            %for now, no coupling of x and v is allowed. v should not appear in delf/delx. for example, x*v^2 not lalowed, atleast for current implementation
               @(x_dot,v_dot)          x_dot+v_dot;...                 %this kind of structure is important, atleasft for the way gradient descent is implemented
               @(x_dot_dot,v_dot_dot)  x_dot_dot+v_dot_dot};            % for now, working only with f^[i] is a function only of ith derivatives of x and v.

g_int_tilde = {@(x,v)                  x-v;
               @(x_dot,v_dot)          -x_dot+v_dot;                 %this kind of structure is important, atleasft for the way gradient descent is implemented
               @(x_dot_dot,v_dot_dot)  x_dot_dot+v_dot_dot};





%---------------------------------------------------------------------------------------------------------------------------------------
%% Initialising required objects
%---------------------------------------------------------------------------------------------------------------------------------------


time=0:dt:T;


tic
disp("Taking Jacobian of the functions...");
jacobian_f_int_tilde_x=jacob(f_int_tilde,["x","x_dot","x_dot_dot"],1);        %do check what this is actually a function of, can make it a .m file later for further clarity, rn all functions would just be functions of the vectore elements w.r.t whose gradient has been sought. 
jacobian_f_int_tilde_v=jacob(f_int_tilde,["v","v_dot","v_dot_dot"],1);        % returns an oxo cell, each having the i-jth elements of jacobian as function handle.

jacobian_g_int_tilde_x=jacob(g_int_tilde,["x","x_dot","x_dot_dot"],1);
jacobian_g_int_tilde_v=jacob(g_int_tilde,["v","v_dot","v_dot_dot"],1);

% jacobian_G_ext_tilde=jacob(G_ext_tilde,["xi_dot","yi_dot"],8); authomation not working, switching to manual mode

% jacobian_G_ext_tilde={@(x,xi,y,yi) -(x-xi)./sqrt((eps+x-xi).^2 + (eps+y-yi).^2), @(x,xi,y,yi) -(y-yi)./sqrt((eps+x-xi).^2 + (eps+y-yi).^2);
%                       @(x,xi,y,yi) 0,@(x,xi,y,yi) 0;
%                       @(x,xi,y,yi) 0,@(x,xi,y,yi) 0}

jacobian_G_ext_tilde={,;
                      ,;
                      ,}
                      



fprintf('Jacobians Ready. Elapsed Time %.2f\n',toc);
disp('Initialising...');






%each of following objects returns a dictionary containing object by their literal names like "RR_tilde"-->{RR_tilde}
action      =init_action(N,d,T);
genprocess  =init_genprocess(o,N,d,T); 
sense       =init_sense(N,o,L,dL,T,gamma_vecdL_y,lambda_vecdL_y);  
genmodel    =init_genmodel(N,o,L,dL,T,gamma_vecdL_x,lambda_vecdL_x,gamma_vecdL_v,lambda_vecdL_v,etas);



%% Change initial configuration. Alternatively change the function init_RR_tilde

RR_tilde=get_from_dict(genprocess,"RR_tilde"); %get the RR_tilde from the dictionary
mu_tilde=get_from_dict(genmodel,"mu_tilde_x"); %get the RR_tilde from the dictionary



% for i=1:Na
    
% end


di=S/(n-1);
for i=1:N
    % RR_tilde{1}(i,:,:)=(2*mod(i,2)-1+i) * ones(1,d,T); 
    if i == 3
        RR_tilde{1}(i,:,1) = [1000, 1000];  % Agent 3 very far
    elseif i == 4
        RR_tilde{1}(i,:,1) = [-1000, -1000]; % Agent 4 very far
    else
        RR_tilde{1}(i,:,1) = [-S/2 + (i-1)*di, 0];  % All other agents aligned along x-axis at y=0
    end
    % RR_tilde{2}(i,:,:)=zeros(1,d,T);    %manually change initial configurations  ||
    % RR_tilde{3}(i,:,:)=zeros(1,d,T);    %      ||
    % mu_tilde{1}{1}(i,:,1)=[10]; 


end
S=100;
genprocess{"RR_tilde"}=RR_tilde;
genmodel{"mu_tilde_x"}=mu_tilde;



%---------------------------------------------------------------------------------------------------------------------------------------
%% simulation loop
%---------------------------------------------------------------------------------------------------------------------------------------


disp('Entering Simulation...');
for t=2:T
    kx=(t>15)*kx0;
    ka=(t>15)*ka0;
    kv=(t>15)*kv0;
if mod(t,10)==0
    fprintf('\rSimulation Left: %.2f %  ', (T-t)*100/T);
end
    
    %each of following objects returns a dictionary containing object by their literal names like "RR_tilde"-->{RR_tilde}
    action      =update_action(genprocess{"RR_tilde"},action,genmodel{"mu_tilde_x"},genmodel{"mu_tilde_v"},g_int_tilde,sense{"Y_ext_tilde"},jacobian_G_ext_tilde,sense{"PI_tilde_y"},ka,dt,t); %to add non internal forces.

                %  update_action(genprocess{"RR_tilde"},genmodel{"action"},genmodel{"mu_tilde_x"},genmodel{"mu_tilde_v"},g_int_tilde,sense{"Y_ext_tilde"},jacobian_G_ext_tilde,PI_tilde_y,t)
    genprocess  =update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"RR_tilde"),...
                                                action,dt,t,S);
    
    sense       =update_sense(sense{"Y_ext_tilde"},G,...
                              genprocess{"RR_tilde"},noise_params,sense{"PI_tilde_y"},...
                              t);  
    
    genmodel    =update_genmodel(genmodel{"vfe"},genmodel{"mu_tilde_x"},genmodel{"mu_tilde_v"},...
                                jacobian_g_int_tilde_x,jacobian_g_int_tilde_v,sense{"PI_tilde_y"},sense{"Y_ext_tilde"},...
                                g_int_tilde,genmodel{"PI_tilde_x"},jacobian_f_int_tilde_x,jacobian_f_int_tilde_v,f_int_tilde,...
                                genmodel{"eta_tilde"},genmodel{"PI_tilde_v"},t,dt,kx,n_kx,kv);
    

    % RR_tilde=get_from_dict(genprocess,"RR_tilde");
    % action(:,:,t)= - RR_tilde{1}(:,:,t-1) ; %To be replaced with genmodel("action")/updateaction
    % genprocess=update_genprocess(get_from_dict(genprocess,"F_tilde"),F,get_from_dict(genprocess,"RR_tilde"),action,dt,t);
    % Y_ext_tilde
end
fprintf('Simulation Complete. Elapsed Time %.2f\n', toc);
disp('PLotting...');

%---------------------------------------------------------------------------------------------------------------------------------------
%% Plotting and analysing
%---------------------------------------------------------------------------------------------------------------------------------------

RR_tilde=genprocess{"RR_tilde"};
vfe=genmodel{"vfe"};
mu=genmodel{"mu_tilde_x"};
Y=sense{"Y_ext_tilde"};


fig = figure('Color', 'w');
mainTitle = sprintf('mu_x, mu_v and action all update dt=%.2f, o=%d, n=%d, N=%d, d=%d, L=%d, dL=%d, T=%d', dt, o, n, N, d, L, dL, T);
paramTitle = sprintf('n_kx=%d, kx=%.2f, kpx=%.2f, kv=%.2f, kpv=%.2f, ka=%.2f', n_kx, kx, kpx, kv, kpv, ka);
tl = tiledlayout(5,1, 'Padding', 'compact', 'TileSpacing', 'compact');
sg = sgtitle(tl, {mainTitle; paramTitle}, 'FontWeight', 'bold', 'FontSize', 12);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

% Define colors
col_agent1 = [0.85 0.33 0.1]; % red
col_agent2 = [0.12 0.45 0.70]; % blue

% 1st subplot: Prediction error
nexttile;
plot(1:T, squeeze(abs(mu{1}{1}(2,1,:) - Y{1}{1}(2,1,:))), 'LineWidth', 1.5, 'Color', col_agent1); hold on;
plot(1:T, squeeze(abs(mu{2}{1}(1,1,:) - Y{2}{1}(1,1,:))), 'LineWidth', 1.5, 'Color', col_agent2);
hold off;
title(sprintf('|Prediction error Vs Time.| \n Agent i sees the agent j'), 'Color', 'k');
xlabel("t");
ylabel('|(g(\mu) - y)|');
legend('Agent 1', 'Agent 2');
grid on;


% 2nd subplot: VFE
nexttile;
plot(1:T, vfe(1,:), 'LineWidth', 1.5, 'Color', col_agent1); hold on;
plot(1:T, vfe(2,:), 'LineWidth', 1.5, 'Color', col_agent2);
hold off;
title("VFE of agent 1 (red) and agent 2 (blue)", 'Color', 'k');
xlabel("t");
ylabel("VFE");
legend('Agent 1', 'Agent 2');
grid on;

% 3rd subplot: Action profile
nexttile;
scatter(1:T, squeeze(action(1,1,:)), 20, col_agent1, 'filled', 'MarkerEdgeColor', col_agent1, 'MarkerFaceAlpha', 0.7); hold on;
scatter(1:T, squeeze(action(2,1,:)), 20, col_agent2, 'filled', 'MarkerEdgeColor', col_agent2, 'MarkerFaceAlpha', 0.7);
hold off;
title("Action profile of agent 1 (red) and 2 (blue), x coordinate", 'Color', 'k');
xlabel("time");
ylabel('Action');
legend('Agent 1', 'Agent 2');
grid on;

% 4th subplot: X coordinates of agent 1 and 2
nexttile;
plot(1:T, squeeze(RR_tilde{1}(1,1,:)), 'LineWidth', 1.5, 'Color', col_agent1); hold on;
plot(1:T, squeeze(RR_tilde{1}(2,1,:)), 'LineWidth', 1.5, 'Color', col_agent2);
hold off;
title('X coordinate of agent 1 (red) and 2 (blue)', 'Color', 'k');
xlabel('time');
ylabel('X position');
legend('Agent 1', 'Agent 2');
grid on;


% 5th subplot: X coordinates of agent 1 and 2
nexttile;
plot(1:T, squeeze(RR_tilde{2}(1,1,:)), 'LineWidth', 1.5, 'Color', col_agent1); hold on;
plot(1:T, squeeze(RR_tilde{2}(2,1,:)), 'LineWidth', 1.5, 'Color', col_agent2);
hold off;
title('Velocity coordinate of agent 1 (red) and 2 (blue)', 'Color', 'k');
xlabel('time');
ylabel('X velocitt');
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

blankfig = figure('Color', 'w');
fig = figure('Color','w'); % black background for visual clarity
% S=max(squeeze(RR_tilde{1}(1,1,:)));

for t = 1:T-10
    clf

    % Plot agent 1 (red) and agent 2 (blue) with specified diameter
    scatter(RR_tilde{1}(1,1,t), RR_tilde{1}(1,2,t), dia, col_agent1, 'filled'); hold on;
    scatter(RR_tilde{1}(2,1,t), RR_tilde{1}(2,2,t), dia, col_agent2, 'filled');
    % scatter(RR_tilde{1}(3,1,t), RR_tilde{1}(2,2,t), dia, 'w', 'filled');
    % scatter(RR_tilde{1}(4,1,t), RR_tilde{1}(2,2,t), dia, 'w', 'filled'); 
    hold off;

    ylabel("Y");
    xlim([-S S]);
    ylim([-S S]);
    set(gca, 'Color', 'k'); % black background for axes
    grid on;
    set(gca, 'GridColor', [1 1 1], 'GridAlpha', 0.5); % white grid lines
    title(sprintf('Time step: %d', t))
    drawnow;

    % Capture the plot as a frame
    frame = getframe(fig);
    writeVideo(outputVideo, frame);
end    
close(outputVideo);
end






ret=RR_tilde;



% end




% update_genprocess(F_tilde,F,RR_tilde,action,dt,t)



















%---------------------------------------------------------------------------------------------------------------------------------------
%% Helper Functions
%---------------------------------------------------------------------------------------------------------------------------------------

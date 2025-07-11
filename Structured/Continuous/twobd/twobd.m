clc
close all;
T=1500;
dt=0.1;
k_mu=0.01*1;
k_a=0.1*1;
alpha=1;
wanna_save=1;
x1=zeros(1,T);
x2=zeros(1,T);

v1=zeros(1,T);
v2=zeros(1,T);

a1=zeros(1,T);
a2=zeros(1,T);

mu_x1=zeros(1,T);
mu_x2=zeros(1,T);

mu_v1=zeros(1,T);
mu_v2=zeros(1,T);

mu_a1=zeros(1,T);
mu_a2=zeros(1,T);

y_x1=zeros(1,T);
y_x2=zeros(1,T);

y_v1=zeros(1,T);
y_v2=zeros(1,T);

grad_mu1=zeros(3,T);
grad_v1=zeros(1,T);

f1=zeros(1,T);
f2=zeros(1,T);


PI_x=init_PI_tilde(3,1,[1],[1]);
PI_y=init_PI_tilde(2,1,[1],[1]);

eta_x=5;

%at t=1
x1(1)=-5;
x2(1)=5;        %putting them on x axis, 10 units apart

v2(1)=0;

mu_x1(1)=eta_x;
mu_x1(2)=eta_x;

y_x1(1)=abs(x2(1)-x1(1));
y_x2(2)=abs(x1(1)-x2(1));

ey_tilde_1=zeros(2,1,T);
ex_tilde_1=zeros(3,1,T);

ey_tilde_2=zeros(2,1,T);
ex_tilde_2=zeros(3,1,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Got gradients from symbolic math %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grad_mu=@(a,eta,ma,mv,mx,yv,yx)[mx.*2.0-yx.*2.0+a.*(mv.*(3.0./2.0)+a.*ma-a.*(eta-mx).*(3.0./2.0))+a.^2.*ma+a.*(mv-a.*(eta-mx)).*(3.0./2.0);mv.*7.0-yv.*4.0+a.*(ma.*2.0+a.*mv.*2.0)+a.*ma.*2.0-a.*(eta-mx).*3.0+a.*(ma+a.*mv).*2.0;ma.*4.0+a.*(mv+a.*ma.*2.0-a.*(eta-mx))+a.*mv.*4.0+a.^2.*ma.*2.0+a.*(mv-a.*(eta-mx))];
grad_Act=@(mv,v1,v2)mv.*4.0+v1.*4.0-v2.*4.0;

%rest things remain 0

% f=@(alpha,mu_x,eta_x) -alpha*(mu_x-eta_x);

for t=2:T
    v1(t)=v1(t-1) + (dt*a1(t-1));
    v2(t)=v2(t-1) + (dt*a2(t-1));

    x1(t)=x1(t-1) + (dt*v1(t));
    x2(t)=x2(t-1) + (dt*v2(t));

    y_x1(t)=abs(x1(t)-x2(t));
    y_v1(t)=v2(t)-v1(t);

    y_x2(t)=abs(x1(t)-x2(t));
    y_v2(t)=v1(t)-v2(t);

    ey_tilde_1(:,:,t)=[y_x1(t) - mu_x1(t-1); 
                       y_v1(t) - mu_v1(t-1)];



    ex_tilde_1(:,:,t)=[mu_v1(t-1) - (-alpha*(mu_x1(t-1)-eta_x));
                       mu_a1(t-1) - (-alpha*(mu_v1(t-1)-0));
                          0       - (-alpha*(mu_a1(t-1)-0))];

    f1(t)= dot(ey_tilde_1(:,:,t),PI_y*ey_tilde_1(:,:,t)) + dot(ex_tilde_1(:,:,t),PI_x*ex_tilde_1(:,:,t));
    



    ey_tilde_2(:,:,t)=[y_x2(t) - mu_x2(t); 
                       y_v2(t) - mu_v2(t)];


    ex_tilde_2(:,:,t)=[mu_v2(t-1) - (-alpha*(mu_x2(t-1)-eta_x));
                       mu_a2(t-1) - (-alpha*(mu_v2(t-1)-0));
                            0     - (-alpha*(mu_a2(t-1)-0))];
    
    f2(t)= dot(ey_tilde_2(:,:,t),PI_y*ey_tilde_2(:,:,t)) + dot(ex_tilde_2(:,:,t),PI_x*ex_tilde_2(:,:,t));



    grad_mu1(:,t)=grad_mu(alpha,eta_x,mu_a1(t-1),mu_v1(t-1),mu_x1(t-1),y_v1(t),y_x1(t));
    
    if t>15
    % mu_x1(t)=mu_x1(t-1) + dt * (mu_v1(t-1)-k_mu*grad(1));
    % mu_v1(t)=mu_v1(t-1) + dt * (mu_a1(t-1)-k_mu*grad(2));
    % mu_a1(t)=mu_a1(t-1) + dt * (0         -k_mu*grad(3));

    mu_x1(t)=mu_x1(t-1) + k_mu * (mu_v1(t-1)-grad_mu1(1,t));
    mu_v1(t)=mu_v1(t-1) + k_mu * (mu_a1(t-1)-grad_mu1(2,t));
    mu_a1(t)=mu_a1(t-1) + k_mu * (0         -grad_mu1(3,t));
    end
    grad1=[0;-1];
    grad_v1(t)=grad_Act(mu_v1(t-1),v1(t),v2(t));
    
    a1(t)=-k_a*grad_v1(t);
    % a2(t)=a2(t-1)+


    
    

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Aesthetic Plotting

% Time vector
time = (1:T) * dt;

% --- Define colors ---
col1 = [0.1 0.6 0.8];   % Agent 1 color (bluish)
col2 = [0.9 0.3 0.3];   % Agent 2 color (reddish)

% --- Define line styles ---
ls1 = '-';
ls2 = '--';

% --- Define parameter title string ---
param_str = sprintf('dt= %.3f,\\alpha = %.3f   k_\\mu = %.3f   k_a = %.3f   \\eta = %.3f',dt, alpha, k_mu, k_a, eta_x);

% --- Create figure ---
figure('Color', 'w');

% === 1. Free Energy ===
subplot(3,2,1)
plot(time, f1, ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, f2, ls2, 'Color', col2, 'LineWidth', 1.5);
title('Variational Free Energy');
xlabel('Time (s)'); ylabel('F');
legend('Agent 1', 'Agent 2');
grid on;

% === 2. Sensory prediction error (position) ===
subplot(3,2,2)
plot(time, squeeze(y_x1(:) - mu_x1(:)), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, squeeze(y_x2(:) - mu_x2(:)), ls2, 'Color', col2, 'LineWidth', 1.5);
title('Prediction Error: Position');
xlabel('Time (s)'); ylabel('\epsilon^y_x');
legend('Agent 1', 'Agent 2');
grid on;

% === 3. Sensory prediction error (velocity) ===
subplot(3,2,3)
plot(time, squeeze(ey_tilde_1(2,1,:)), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, squeeze(ey_tilde_2(2,1,:)), ls2, 'Color', col2, 'LineWidth', 1.5);
title('Prediction Error: Velocity');
xlabel('Time (s)'); ylabel('\epsilon^y_v');
legend('Agent 1', 'Agent 2');
grid on;

% === 4. Hidden state: x1 ===
subplot(3,2,4)
plot(time, squeeze(x1), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
% plot(time, squeeze(x2), ls2, 'Color', col2, 'LineWidth', 1.5); % Uncomment for Agent 2
title('Hidden State: x_1');
xlabel('Time (s)'); ylabel('x_1');
legend('Agent 1');
grid on;

% === 5. Action gradient (acceleration proxy) ===
subplot(3,2,5)
scatter(time, grad_v1, 15, col1, 'filled'); hold on;
% scatter(time, grad_v2, 15, col2, 'filled'); % Uncomment for Agent 2
title('Action Gradient (acceleration)');
xlabel('Time (s)'); ylabel('d(a)/dt \approx \epsilon^x_2');
legend('Agent 1');
grid on;

% === 6. (Optional) Hidden prediction error: jerk ===
% subplot(3,2,6)
% plot(time, squeeze(ex_tilde_1(3,1,:)), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
% plot(time, squeeze(ex_tilde_2(3,1,:)), ls2, 'Color', col2, 'LineWidth', 1.5);
% title('Hidden State Error: Jerk');
% xlabel('Time (s)'); ylabel('\epsilon^x_3');
% legend('Agent 1', 'Agent 2');
% grid on;

% --- Title for all subplots ---
sgtitle({'Two-body Active Inference Simulation', ...
         ['\fontsize{10}', param_str]}, ...
         'FontWeight', 'bold', 'FontSize', 14);











if wanna_save
output_dir ="D:\Projects\Summer 2025\Surprise Minimisation\outputs\structured_testing"; % replace with your desired path

now_dt = datetime('now');
now_dt.TimeZone = 'local';
now=string(now_dt, 'd-MMM-y HH:mm:ss');
now=split(now,":");
now=join(now,"_");
filename = fullfile(output_dir, sprintf('twobd_%s', ...
now));
% Define colors
col_agent1 = [0.85 0.33 0.1]; % red
col_agent2 = [0.12 0.45 0.70]; % blue
S=10;
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
    scatter(x1(t), 0, 100, col_agent1, 'filled'); hold on;
    scatter(x2(t), 0, 100, col_agent1, 'filled');
    % scatter(RR_tilde{1}(3,1,t), RR_tilde{1}(2,2,t), dia, 'w', 'filled');
    % scatter(RR_tilde{1}(4,1,t), RR_tilde{1}(2,2,t), dia, 'w', 'filled'); 
    hold off;
    xlabel("X")
    ylabel("Y");
    xlim([-S S]);
    ylim([-10 10]);
    set(gca, 'Color', 'k'); % black background for axes
    grid on;
    set(gca, 'GridColor', [1 1 1], 'GridAlpha', 0.5); % white grid lines
    title(sprintf('t= %d, dt= %.2f, alpha= %.2f, kmu= %.2f, ka=%.2f, eta=%.2f',t,dt,alpha,k_mu,k_a,eta_x))
    drawnow;

    % Capture the plot as a frame
    frame = getframe(fig);
    writeVideo(outputVideo, frame);
end    
close(outputVideo);
end




















function ret=init_PI_tilde(o,L,gamma_vecdL,lambda_vecdL)
    %return a cell(1,dl), each cell containing olxol matrix. 

    %gamma and lambda _vecdl are 1xdl vectors
    dL=length(gamma_vecdL); %by default dL=1
    PI=cell(1,dL);

    parfor j=1:dL
        Spatial_term_j=inv(gamma_vecdL(j)*eye(L));
        Temporal_term_j=spm_DEM_R(o,lambda_vecdL(j));
        % PI{j}=kron(Spatial_term_j,Temporal_term_j);
        PI{j}=kron(Temporal_term_j,Spatial_term_j); %kept like this intentionally, w.r.t how I just flatten out mu_tilde later
    end

    ret=PI{1};
end









































clc
close all;
T=15;
dt=0.1;
to_mu=0;
to_act=0;
k_mu=0.01*1;
k_a=0.1*1;
N=2;
L=1;
alpha=1;
wanna_save=1;
v0=1;
T=100;
L=1;

%agents senses it's own angle theta, v0 is fixed


pos=zeros(2,2,T);
vel=zeros(2,2,T);
acc=zeros(2,2,T);


mu_theta=zeros(2,T);
mu_theta_dot=zeros(2,T);
mu_theta_dot_dot=zeros(2,T);

theta=zeros(2,T);
theta_dot=zeros(2,T);

grad_mu=zeros(2,3,T);
grad_vel=zeros(2,2,T);

fe=zeros(2,T);

PI_x=init_PI_tilde(3,L,[1],[1]);
PI_y=init_PI_tilde(2,L,[1],[1]);

eta_x=0;

pos(1,1,1)=5;
pos(1,2,1)=1;
pos(2,1,1)=-5;
pos(2,2,1)=-1;

vel(:,:,1)=randn(2,2);

mu_theta(:,1)=eta_x;

theta(1,1)=atan2d(vel(1,2,1),vel(1,1,1));
theta(2,1)=atan2d(vel(2,2,1),vel(2,1,1));

ey_tilde_1=zeros(2,1,T);
ex_tilde_1=zeros(3,1,T);

ey_tilde_2=zeros(2,1,T);
ex_tilde_2=zeros(3,1,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Got gradients from symbolic math %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GM=@(alpha,eta_x,mu_theta,mu_theta_dot,theta,theta_dot)[mu_theta.*6.0-theta.*2.0+alpha.*(mu_theta_dot+alpha.*mu_theta.*2.0-alpha.*(eta_x-mu_theta))-alpha.*(eta_x-mu_theta_dot).*4.0+alpha.*(mu_theta_dot.*(3.0./2.0)+alpha.*mu_theta-alpha.*(eta_x-mu_theta).*(3.0./2.0))+alpha.^2.*mu_theta.*3.0+alpha.*(mu_theta_dot-alpha.*(eta_x-mu_theta)).*(5.0./2.0);mu_theta_dot.*7.0-theta_dot.*4.0+alpha.*mu_theta.*2.0+alpha.*(mu_theta.*2.0-alpha.*(eta_x-mu_theta_dot).*2.0)-alpha.*(eta_x-mu_theta).*3.0+alpha.*(mu_theta-alpha.*(eta_x-mu_theta_dot)).*2.0;mu_theta.*6.0-theta.*2.0+alpha.*(mu_theta_dot+alpha.*mu_theta.*2.0-alpha.*(eta_x-mu_theta))-alpha.*(eta_x-mu_theta_dot).*4.0+alpha.*(mu_theta_dot.*(3.0./2.0)+alpha.*mu_theta-alpha.*(eta_x-mu_theta).*(3.0./2.0))+alpha.^2.*mu_theta.*3.0+alpha.*(mu_theta_dot-alpha.*(eta_x-mu_theta)).*(5.0./2.0)]
% GA=

for t=2:T
    vel(:,:,t)=vel(:,:,t-1) +(dt*acc(:,:,t-1));

    pos(:,:,t)=pos(:,:,t-1) +(dt*vel(:,:,t));

    theta(1,t)=atan2d(vel(1,2,t),vel(1,1,t));
    theta(2,t)=atan2d(vel(2,2,t),vel(2,1,t));

    theta_dot(1,t) = (vel(1,1,t) * acc(1,2,t) - vel(1,2,t) * acc(1,1,t)) / (eps + vel(1,1,t)^2 + vel(1,2,t)^2);
    theta_dot(2,t) = (vel(2,1,t) * acc(2,2,t) - vel(2,2,t) * acc(2,1,t)) / (eps + vel(2,1,t)^2 + vel(2,2,t)^2);

    ey_tilde_1(:,:,t)=[theta(1,t)    - mu_theta(1,t-1);
                       theta_dot(1,t)- mu_theta_dot(1,t-1)];
    
    ex_tilde_1(:,:,t)=[mu_theta_dot(1,t-1)     - (-alpha*(mu_theta(1,t-1)-eta_x));
                       mu_theta_dot_dot(1,t-1) - (-alpha*(mu_theta_dot(1,t-1)-eta_x));
                                               - (-alpha*(mu_theta_dot_dot(1,t-1)-0))];

    fe(1,t)=dot(ey_tilde_1(:,:,t),PI_y*ey_tilde_1(:,:,t)) + dot(ex_tilde_1(:,:,t),PI_x*ex_tilde_1(:,:,t));



    ey_tilde_2(:,:,t)=[theta(2,t)    - mu_theta(2,t-1);
                       theta_dot(2,t)- mu_theta_dot(2,t-1)];
    
    ex_tilde_2(:,:,t)=[mu_theta_dot(2,t-1)     - (-alpha*(mu_theta(2,t-1)-eta_x));
                       mu_theta_dot_dot(2,t-1) - (-alpha*(mu_theta_dot(2,t-1)-eta_x));
                                               - (-alpha*(mu_theta_dot_dot(2,t-1)-0))];

    fe(2,t)=dot(ey_tilde_2(:,:,t),PI_y*ey_tilde_2(:,:,t)) + dot(ex_tilde_2(:,:,t),PI_x*ex_tilde_2(:,:,t));

    grad_mu(1,:,t)=GM(alpha,eta_x,mu_theta(1,t-1),mu_theta_dot(1,t-1),theta(1,t),theta_dot(1,t))
    grad_mu(2,:,t)=GM(alpha,eta_x,mu_theta(2,t-1),mu_theta_dot(2,t-1),theta(2,t),theta_dot(2,t))

    mu_theta(1,t)        =(t>to_mu).* (mu_theta(1,t-1)         + k_mu *(mu_theta_dot(t-1)    -grad_mu(1,1,t)));
    mu_theta_dot(1,t)    =(t>to_mu).* (mu_theta_dot(1,t-1)     + k_mu *(mu_theta_dot_dot(t-1)-grad_mu(1,2,t)));
    mu_theta_dot_dot(1,t)=(t>to_mu).* (mu_theta_dot_dot(1,t-1) + k_mu *(0                    -grad_mu(1,3,t)));

end




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Aesthetic Plotting

% Time vector
time = (1:T);

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
plot(time, fe(1,:), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, fe(1,:), ls2, 'Color', col2, 'LineWidth', 1.5);
title('Variational Free Energy');
xlabel('Time (s)'); ylabel('F');
legend('Agent 1', 'Agent 2');
grid on;

% === 2. Sensory prediction error (position) ===
subplot(3,2,2)
plot(time, squeeze(theta(1,:) - mu_theta(1,:)), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, squeeze(theta(2,:) - mu_theta(2,:)), ls2, 'Color', col2, 'LineWidth', 1.5);
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
plot(time, squeeze(pos(1,1,:)), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, squeeze(pos(2,1,:)), ls2, 'Color', col2, 'LineWidth', 1.5); % Uncomment for Agent 2
title('x_1');
xlabel('Time (s)'); ylabel('x_1');
legend('Agent 1');
grid on;

% === 5. Action gradient (acceleration proxy) ===
subplot(3,2,5)
scatter(time, squeeze(grad_mu(1,1,:)), 15, col1, 'filled'); hold on;
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













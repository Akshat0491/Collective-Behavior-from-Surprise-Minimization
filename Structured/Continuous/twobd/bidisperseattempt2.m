close all;

N=2;
T=1000;
dt=0.1;
k_mu=0.1*1;
k_a=1*1;
alph=1;
isnormed=1;
wanna_save=1;
to_mu=0;
v0=1;

v0x=1;
v0y=0;

range_of_interaction=10;

nsims=1;
thetass=zeros(nsims,T);


% parfor n=1:nsims


pos=zeros(N,2,T);
vel=zeros(N,2,T);
acc=zeros(N,2,T);

dummy_vector=[v0x,v0y];

y_phi=zeros(N,T);
% y_phi_dot=zeros(N,T);

mu_phi=zeros(N,T);
mu_phi_dot=zeros(N,T);

noise_params=[1,1,1];
PI_tilde_x=diag(noise_params(1:2));
PI_tilde_y=noise_params(3);

fe=zeros(2,T);

eta0=1;

S=20;

% rng(42);  % for reproducibility


dists=zeros(1,T);

pos(1,1,1)=5;
pos(1,2,1)=5;
pos(2,1,1)=-5;
pos(2,2,1)=0;

% vel(:,:,1)=randn(N,2);
% vel(:,:,1)=ones(N,2);
vel(1,:,1)=[1,0];
vel(2,:,1)=[-1,0];
vel=normalise_vel(vel,v0,1);

mu_phi(:,1)=eta0;


ey_tilde=zeros(2,N,T);
ex_tilde=zeros(2,N,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Gradients
grad_mu=@(alph,eta_x,mu_phi,mu_phi_dot,y_phi)[mu_phi.*2.0-y_phi.*2.0+alph.*(mu_phi_dot-alph.*(eta_x-mu_phi)).*2.0;mu_phi_dot.*2.0-alph.*(eta_x-mu_phi).*2.0+alph.^2.*mu_phi_dot.*2.0];
gradmuu=cell(N,1);
grad_act=@(mu_phi,v0x,v0y,vx,vy)[v0x.*(-mu_phi+v0x.*vx+v0y.*vy).*2.0;v0y.*(-mu_phi+v0x.*vx+v0y.*vy).*2.0]
gradactt=cell(N,1);

% ey_tilde1=zeros(2,T);

% ex_tilde1=zeros(2,N,T);

for i=1:N
    eta=((-1)^(2-i))*eta0;
    gradmuu{i}=zeros(2,T);
    gradactt{i}=zeros(2,T);
    % gradactt{i}=zeros(1,T)
    y_phi(i,1)=dot(vel(i,:,1),dummy_vector);
    % y_phi_dot(i,1)=

    ey_tilde(:,i,1)=[y_phi(i,1)-mu_phi(i,1)];
    % ;y_phi_dot(i,1)-mu_phi_dot(i,1)]

    ex_tilde(:,i,1)=[mu_phi_dot(i,1)-(-alph*(mu_phi(i,1)-1));
                     -(-alph*mu_phi_dot(i,1))];

    % fe(i,1)=dot(ey_tilde(:,i,1),PI_tilde_y*ey_tilde(:,i,1))+dot(ex_tilde(:,i,1),PI_tilde_x*ex_tilde(:,i,1));
end



for t=2:T
    for i=1:N
            vel(i,:,t)=vel(i,:,t-1) + dt*acc(i,:,t-1);
    vel=normalise_vel(vel,v0,t);
    pos(i,:,t)=pos(i,:,t-1) + dt*vel(i,:,t);


            
    pos(i,:,t)=pos(i,:,t);

    % pos(:,:,t) = pos(:,:,t-1) + dt * vel(:,:,t);

% Apply periodic boundary around [-S, S] for both X and Y
pos(i,:,t) = mod(pos(i,:,t) + S, 2*S) - S;

    
        eta=((-1)^(2-i))*eta0;
        y_phi(i,t)=dot(vel(i,:,t),dummy_vector);
        ey_tilde(:,i,t)=[y_phi(i,t)-mu_phi(i,t-1)];
    % ;y_phi_dot(i,1)-mu_phi_dot(i,1)]

        ex_tilde(:,i,t)=[mu_phi_dot(i,t-1)-(-alph*(mu_phi(i,t-1)-1));
                         -(-alph*mu_phi_dot(i,t-1))];


        fe(i,t)=dot(ey_tilde(:,i,t),PI_tilde_y*ey_tilde(:,i,t))+dot(ex_tilde(:,i,t),PI_tilde_x*ex_tilde(:,i,t));


        gradmuu{i}(:,t)=grad_mu(alph,eta,mu_phi(i,t-1),mu_phi_dot(i,t-1),y_phi(i,t));


        mu_phi(i,t)    =(t>to_mu).* (mu_phi(i,t-1)     + k_mu * (mu_phi_dot(i,t-1)-gradmuu{i}(1,t)));
        mu_phi_dot(i,t)=(t>to_mu).* (mu_phi_dot(i,t-1) + k_mu * (0-gradmuu{i}(2,t)));

        % acc(i,:,t)=-(1).*(1*k_a/dt).*(y_phi(i,t)-mu_phi(i,t)).*(dummy_vector);

        

        acc(i,:,t)=(-10*k_a*grad_act(mu_phi(i,t),v0x,v0y,vel(i,1,t),vel(i,2,t)));
    end


end

% thetass(n,:)=squeeze(atan2d(vel(1,2,:),vel(1,1,:)));
% n

% end



























































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
param_str = sprintf('dt= %.3f,\\alph = %.3f   k_\\mu = %.3f   k_a = %.3f   \\eta = %.3f isnormed= %d ',dt, alph, k_mu, k_a, eta,isnormed);

% --- Create figure ---
figure('Color', 'w');

% === 1. Free Energy ===
subplot(3,2,1)
plot(time, fe(1,:), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, fe(2,:), ls2, 'Color', col2, 'LineWidth', 1.5);
title('Variational Free Energy');
xlabel('Time (s)'); ylabel('F');
legend('Agent 1', 'Agent 2');
grid on;

% === 2. Sensory prediction error (position) ===
subplot(3,2,2)
plot(time, squeeze(y_phi(1,:) - mu_phi(1,:)), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, squeeze(y_phi(2,:) - mu_phi(2,:)), ls2, 'Color', col2, 'LineWidth', 1.5);
title('Prediction Error: cos(y_\phi)-cos(mu_\phi)');
xlabel('Time (s)'); ylabel('\epsilon_y');
legend('Agent 1', 'Agent 2');
grid on;

% === 3. Sensations ===
subplot(3,2,3)
plot(time, squeeze(y_phi(1,:)), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
% plot(time, squeeze(ey_tilde_2(2,1,:)), ls2, 'Color', col2, 'LineWidth', 1.5);
title('cos(\phi)');
xlabel('Time (s)'); ylabel('cos(\phi)');
legend('Agent 1');
grid on;

% === 4. Hidden state: x1 ===
subplot(3,2,4)
plot(time, squeeze(atan2d(vel(1,2,:),vel(1,1,:))), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
plot(time, squeeze(atan2d(vel(2,2,:),vel(2,1,:))), ls2, 'Color', col2, 'LineWidth', 1.5); % Uncomment for Agent 2
title('x_1');
xlabel('Time (s)'); ylabel('x_1');
legend('Agent 1');
grid on;

% === 5. Action gradient (acceleration proxy) ===
subplot(3,2,5)
scatter(time, squeeze(acc(1,1,:)), 15, col1, 'filled'); hold on;
% scatter(time, grad_v2, 15, col2, 'filled'); % Uncomment for Agent 2
title('Action Gradient (acceleration)');
xlabel('Time (s)'); ylabel('d(a)/dt \approx \epsilon^x_2');
legend('Agent 1');
grid on;

% % === 6. (Optional) Hidden prediction error: jerk ===
% % subplot(3,2,6)
% % plot(time, squeeze(ex_tilde_1(3,1,:)), ls1, 'Color', col1, 'LineWidth', 1.5); hold on;
% % plot(time, squeeze(ex_tilde_2(3,1,:)), ls2, 'Color', col2, 'LineWidth', 1.5);
% % title('Hidden State Error: Jerk');
% % xlabel('Time (s)'); ylabel('\epsilon^x_3');
% % legend('Agent 1', 'Agent 2');
% % grid on;

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

outputVideo = VideoWriter(filename,'MPEG-4');
outputVideo.FrameRate = 70; % Adjust as needed
% outputVideo.FileFormat='mp4';
open(outputVideo);

blankfig = figure('Color', 'w');
fig = figure('Color','w'); % black background for visual clarity
% S=max(squeeze(RR_tilde{1}(1,1,:)));

for t = 1:T-10
    clf

    % Plot agent 1 (red) and agent 2 (blue) with specified diameter
    scatter(pos(1,1,t), pos(1,2,t), 100, col_agent1, 'filled'); hold on;
    quiver(pos(1,1,t), pos(1,2,t),vel(1,1,t), vel(1,2,t),'filled')
    scatter(pos(2,1,t), pos(2,2,t), 100, col_agent2, 'filled');
    quiver(pos(2,1,t), pos(2,2,t),vel(2,1,t), vel(2,2,t),'filled')
    % scatter(RR_tilde{1}(3,1,t), RR_tilde{1}(2,2,t), dia, 'w', 'filled');
    % scatter(RR_tilde{1}(4,1,t), RR_tilde{1}(2,2,t), dia, 'w', 'filled'); 
    hold off;
    xlabel("X")
    ylabel("Y");
    xlim([-(30) (30)]);
    ylim([-30 30]);
    % ylim padded; xlim padded;

    set(gca, 'Color', 'w'); % black background for axes
    grid on;
    set(gca, 'GridColor', [1 1 1], 'GridAlpha', 0.5); % white grid lines
    title(sprintf('t= %d, dt= %.2f, alpha= %.2f, kmu= %.2f, ka=%.2f, eta=%.2f, isnormed= %d',t,dt,alph,k_mu,k_a,eta,isnormed))
    drawnow;

    % Capture the plot as a frame
    frame = getframe(fig);
    % writeVideo(outputVideo, frame);
end    
close(outputVideo);
end





% thetass(n,:)=squeeze(atan2d(vel(1,2,:),vel(1,1,:)));
% end



% plot(1:T,thetass)








function ret=get_angle(i,vel,t)
    ret=atan2(vel(i,2,t),vel(i,1,t));
end

function ret=normalise_vel(vel,v0,t)
    temp=vel(:,:,t); %Nx2
    nor=sqrt(sum(temp.^2,2));
    vel(:,:,t)=v0.*temp./nor;
    ret=vel;
end

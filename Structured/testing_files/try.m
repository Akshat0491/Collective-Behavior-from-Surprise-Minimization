clc;
clear;
close all;

%% Parameters
m = 1;        % mass
k = 1;        % spring constant
kc = 0.2;     % coupling spring constant
t_max = 50;   % simulation time
dt = 0.01;    % time step

%% Time vector
t = 0:dt:t_max;
N = length(t);

%% Preallocate
x1 = zeros(1, N); v1 = zeros(1, N);
x2 = zeros(1, N); v2 = zeros(1, N);

%% Initial conditions
x1(1) = 1;   % initial displacement oscillator 1
v1(1) = 0;   % initial velocity oscillator 1
x2(1) = 0;   % initial displacement oscillator 2
v2(1) = 0;   % initial velocity oscillator 2

%% Integration (Euler method)
for i = 1:N-1
    % Accelerations
    a1 = -(k*x1(i) + kc*(x1(i) - x2(i))) / m;
    a2 = -(k*x2(i) + kc*(x2(i) - x1(i))) / m;

    % Euler updates
    v1(i+1) = v1(i) + a1*dt;
    x1(i+1) = x1(i) + v1(i)*dt;

    v2(i+1) = v2(i) + a2*dt;
    x2(i+1) = x2(i) + v2(i)*dt;
end

%% Plot Phase Space Diagrams

figure;
subplot(2,2,1)
plot(x1, v1, 'b')
xlabel('x_1'); ylabel('v_1');
title('Phase Space: Oscillator 1');
grid on

subplot(2,2,2)
plot(x2, v2, 'r')
xlabel('x_2'); ylabel('v_2');
title('Phase Space: Oscillator 2');
grid on

subplot(2,2,[3 4])
plot(x1, x2, 'k')
xlabel('x_1'); ylabel('x_2');
title('x_1 vs x_2 (Coupling Trajectory)');
grid on


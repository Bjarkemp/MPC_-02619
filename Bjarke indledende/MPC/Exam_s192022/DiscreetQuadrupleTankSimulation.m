clc; clear all; close all; clearvars;

% Main simulation script without control
t0 = 0;
tf = 1200; % 20 minutes in seconds
Ts = 10; % Sampling time [s]
t = t0:Ts:tf; % Sampling instants [s]
num_steps = length(t);

% Initial parameters
m10 = 0.0; m20 = 0.0; m30 = 0.0; m40 = 0.0;
x0 = [m10; m20; m30; m40];
F1 = 300; F2 = 300;
u = [repmat(F1,1,num_steps); repmat(F2,1,num_steps)];
umin = [0; 0];
umax = [400; 1000];

% Parameters
p = [1.2272; 1.2272; 1.2272; 1.2272; 380.1327; 380.1327; 380.1327; 380.1327; 981; 0.45; 0.40; 1];

% Pre-allocate for performance
num_steps = length(t);
X = zeros(num_steps, 4); % System states
T = zeros(num_steps, 1); % Time vector
U = zeros(num_steps, 2); % Control inputs (F1 and F2)

y = zeros(4, num_steps); % Measured states (with noise)
z = zeros(2, num_steps); % Simplified measured outputs (m1 and m2)

% Process Noise
Q = [20^2 0; 0 40^2]; % Process noise covariance matrix
Lq = chol(Q,'lower');
w = Lq * randn(2, num_steps); % Process noise for each time step

% Measurement Noise
R = eye(4); % Measurement noise covariance matrix
Lr = chol(R,'lower');
v = Lr * randn(4, num_steps); % Measurement noise for each time step

% Initial state
x = x0;

for k = 1:num_steps-1
    % Measurement with noise
    y(:,k) = x + v(:,k); % Measured levels with noise at time step k
    z(:,k) = [x(1), x(2)]; % Simplified measurement for first two tanks
    
    % Control input (no controller in this case, keeping constant inputs F1 and F2)
    %u = [F1; F2]; % Using constant inputs in this simplified case
    
    % Simulate process from t(k) to t(k+1)
    [T_temp, X_temp] = ode15s(@(t,x) QuadrupleTankProcess(t, x, u + w(:,k), p), [t(k) t(k+1)], x);
    
    % Update state
    x = X_temp(end, :)'; % Final state after time step
    T(k+1) = t(k+1); % Update time record
    X(k+1, :) = x'; % Store state history   
end

% Plot results
figure;
% Plot tank liquid mass (system response)
subplot(2,1,1);
plot(T, X);
xlabel('Time (s)');
ylabel('Tank liquid mass (g)');
legend('m1', 'm2', 'm3', 'm4');
title('Quadruple Tank System Response (Without Control)');

% Plot manipulated variables (control inputs)
subplot(2,1,2);
stairs(T, (u+w)');
xlabel('Time (s)');
ylabel('Control Input (cm^3/s)');
legend('F1', 'F2');
title('Control Inputs over Time');



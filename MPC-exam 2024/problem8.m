clc, clear, close all
addpath("Functions");

% Problem 5
% ----------------------------------------------------------
% Parameters
% ----------------------------------------------------------
a1 = 1.2272;                % [cm2] Area of outlet pipe (Tank 1)
a2 = 1.2272;                % [cm2] Area of outlet pipe (Tank 2)
a3 = 1.2272;                % [cm2] Area of outlet pipe (Tank 3)
a4 = 1.2272;                % [cm2] Area of outlet pipe (Tank 4)
A1 = 380.1327;              % [cm2] Cross sectional area (Tank 1)
A2 = 380.1327;              % [cm2] Cross sectional area (Tank 2)
A3 = 380.1327;              % [cm2] Cross sectional area (Tank 3)
A4 = 380.1327;              % [cm2] Cross sectional area (Tank 4)
g = 981;                    % [cm/s2] Gravity
gamma1 = 0.6;               % Flow distribution constant (valve 1)
gamma2 = 0.7;               % Flow distribution constant (valve 2)
rho = 1.0;                  % [g/cm^3] Density of water

% Parameter-vector
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho]; 
at = p(1:4);                % [cm2] Area of outlet pipe
At = p(5:8);                % [cm2] Cross sectional area

% -----------------------------------------------------------
% Simulation scenario
% -----------------------------------------------------------
t0 = 0.0;                   % [s] Initial time
tf= 20*60;                  % [s] End time
dt = 20;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
Ph = 5;                     % Prediction horizon
m10 = 17612.0123864868;                    % [g] Liquid mass in tank 1 at time t0
m20 = 29640.6694933624;                    % [g] Liquid mass in tank 2 at time t0
m30 = 4644.21948249842;                    % [g] Liquid mass in tank 3 at time t0
m40 = 9378.49308238599;                    % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 = 100;
F4_0 = 150;
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];           % [cm3/s] Disturbance variables at t0
d = d0.*ones(2, length(t));
u = u0.*ones(2, length(t));
[y0] = sensor_wo_noise(x0', At, rho);

%%  -------------------- 8.1 MPC function ------------------------------

%linearization
% Steady State
xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);
ys = sensor_wo_noise(xs,at,rho);
zs = sensor_wo_noise(xs,at,rho);

%Stochastic Brownian
R = [1^2 0 0 0; 0 1^2 0 0; 0 0 0.5^2 0; 0 0 0 0.5^2];     % Covariance for measurement noise
sigma = [2^2 0; 0 2^2]; 
[A, B, C, G, Gw] = linearized_models(xs, At, at, rho, gamma1, gamma2, g, ...
                                                                                'brownian', d0, sigma, R);
%ZOH Discretization of Linear System
%Stochastic Brownian
[Ad,Bd,Gd]=c2dzoh(A,B,G,dt);
D = zeros(2,4);
sys = ss(Ad,[Bd,Gd],C(1:2,:),D);

A = sys.A;
B = sys.B;
C = sys.C;

% %Markov parameters
% [x11_brownian, x12_brownian, x21_brownian, x22_brownian] = MarkovPara(Ad,Bd,C,D,N);

%generate constants (this is moved inside "UnconstrainedMPCDesign")
% phi = generate_phi(A, C, Ph);
% Gamma = generate_Gamma(A, B, C, Ph);

% Define MPC parameters
Q = 100 * eye(size(C, 1));  % Weight on output tracking
S = 0.1 * eye(size(B, 2)); % Weight on control effort
% N is Prediction horizon

% Design MPC
MPC_sys = UnconstrainedMPCDesign(A, B, C, Q, S, Ph);

% Kalman filter parameters
R_hat = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4;     % Covariance for measurement noise
Q_hat = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;           % Covariance for process noise

u_vec = zeros(1, tf / dt); % Initial control input
R = repmat(ones(1, Ph), tf / dt, 1); % Constant reference trajectory

% Generate noisy measurements
y_meas = C * x0 + 0.1 * randn(1, tf / dt); % Example noisy measurement

% Pack inputs
R_full = repmat(ones(size(C,1), 1), Ph + tf/dt, 1);
inputs = {x0, u_vec, R, y_meas};

% Simulation
[y, u, x_hat] = MPC_Sim_Unconstrained(A, B, C, MPC_sys, Q_hat, R_hat, tf, dt, inputs, Ph)

% Plot results
t = (0:dt:tf-dt);
figure;
subplot(3, 1, 1);
plot(t, y);
title('Output Trajectory');
xlabel('Time (s)');
ylabel('Output');

subplot(3, 1, 2);
plot(t, u);
title('Control Input');
xlabel('Time (s)');
ylabel('Input');

subplot(3, 1, 3);
plot(t, x_hat(1, 1:end-1));
title('State Estimate (x1)');
xlabel('Time (s)');
ylabel('State');
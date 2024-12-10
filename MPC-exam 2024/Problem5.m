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
dt = 1;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
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

%%  -------------------- 5.1 Linearization ------------------------------
% Steady State
xs0 = [5000; 5000; 5000; 5000]; % [g] Initial guess on xs
xs = fsolve(@FourTankSystemWrap,xs0,[],u0,d0,p);
ys = sensor_wo_noise(xs,at,rho);
zs = sensor_wo_noise(xs,at,rho);

%Deterministic
[A_det, B_det, C_det, Cz_det, G_det, Gw_det] = system_matrices(xs, At, at, rho, gamma1, gamma2, g);

%Stochastic
R = [1^2 0 0 0; 0 1^2 0 0; 0 0 0.5^2 0; 0 0 0 0.5^2];     % Covariance for measurement noise
[A_stoch, B_stoch, C_stoch, G_stoch, Gw_stoch] = linearized_models(xs, At, at, rho, gamma1, gamma2, g, ...
                                                                  'stochastic', [], [], R);
%Stochastic Brownian
sigma = [2^2 0; 0 2^2]; 
[A_brownian, B_brownian, C_brownian, G_brownian, Gw_brownian] = linearized_models(xs, At, at, rho, gamma1, gamma2, g, ...
                                                                                'brownian', d0, sigma, R);
%%  -------------------- 5.2 Transfer functions ------------------------------
% Define D matrix (assumeing zeros)
D = zeros(size(C_det, 1), size(B_det, 2));

%Deterministic
G_tf_det = compute_transfer_functions(A_det, B_det, C_det, D)

%Stochastic
G_tf_stoch = compute_transfer_functions(A_stoch, B_stoch, C_stoch, D)

%Stochastic brownian
G_tf_brownian = compute_transfer_functions(A_brownian, B_brownian, C_brownian, D)

%%  -------------------- 5.3 Comparison to problem 4 ------------------------------

%take values from problem 4 (K1, K2, tau_1, tau_2) compare to following

[gains_det, time_constants_det] = analyze_transfer_functions(G_tf_det)
[gains_stoch, time_constants_stoch] = analyze_transfer_functions(G_tf_stoch)
[gains_brownian, time_constants_brownian] = analyze_transfer_functions(G_tf_brownian)

%%  -------------------- 5.4 discrete-time state space models ------------------------------
%ZOH Discretization of Linear System
%Deterministic
[Ad_det,Bd_det,Gd_det]=c2dzoh(A_det,B_det,G_det,dt);

%Stochastic 
[Ad_stoch,Bd_stoch,Gd_stoch]=c2dzoh(A_stoch,B_stoch,G_stoch,dt);

%Stochastic Brownian
[Ad_brownian,Bd_brownian,Gd_brownian]=c2dzoh(A_brownian,B_brownian,G_brownian,dt);

%%  -------------------- 5.5 Markov parameters ------------------------------

%sammenlign med markov parametrene fra 4 ved at lave plots

[x11_det, x12_det , x21_det, x22_det] = MarkovPara(Ad_det,Bd_det,C_det,D,N);
[x11_stoch, x12_stoch, x21_stoch, x22_stoch] = MarkovPara(Ad_stoch,Bd_stoch,C_stoch,D,N);
[x11_brownian, x12_brownian, x21_brownian, x22_brownian] = MarkovPara(Ad_brownian,Bd_brownian,C_brownian,D,N);

% figure;
% hold on;
% plot(t/60, x11_det, '-.k', 'LineWidth', 1.5);
% plot(t/60, x12_det, '--c', 'LineWidth', 1.5);
% plot(t/60, x21_det, ':m', 'LineWidth', 1.5);
% plot(t/60, x22_det, '-g', 'LineWidth', 1.5);
% hold off;
% legend('u_1 \rightarrow y_1', 'u_1 \rightarrow y_2', 'u_2 \rightarrow y_1', 'u_2 \rightarrow y_2', ...
%     'Location', 'northwest');
% xlabel('Time [min]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Height [cm]', 'FontSize', 12, 'FontWeight', 'bold'); 
% title('Response of System Outputs to Inputs', 'FontSize', 14, 'FontWeight', 'bold');
% 
% figure;
% hold on;
% plot(t/60, x11_stoch, '-.k', 'LineWidth', 1.5);
% plot(t/60, x12_stoch, '--c', 'LineWidth', 1.5);
% plot(t/60, x21_stoch, ':m', 'LineWidth', 1.5);
% plot(t/60, x22_stoch, '-g', 'LineWidth', 1.5);
% hold off;
% legend('u_1 \rightarrow y_1', 'u_1 \rightarrow y_2', 'u_2 \rightarrow y_1', 'u_2 \rightarrow y_2', ...
%     'Location', 'northwest');
% xlabel('Time [min]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Height [cm]', 'FontSize', 12, 'FontWeight', 'bold'); 
% title('Response of System Outputs to Inputs', 'FontSize', 14, 'FontWeight', 'bold');

figure;
hold on;
plot(t/60, x11_brownian, '-.k', 'LineWidth', 1.5);
plot(t/60, x12_brownian, '--c', 'LineWidth', 1.5);
plot(t/60, x21_brownian, ':m', 'LineWidth', 1.5);
plot(t/60, x22_brownian, '-g', 'LineWidth', 1.5);
hold off;
legend('u_1 \rightarrow y_1', 'u_1 \rightarrow y_2', 'u_2 \rightarrow y_1', 'u_2 \rightarrow y_2', ...
    'Location', 'northwest');
xlabel('Time [min]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Height [cm]', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Response of System Outputs to Inputs', 'FontSize', 14, 'FontWeight', 'bold');

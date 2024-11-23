clc, clear, close all
addpath("Functions");
%% Problem 2
% ----------------------------------------------------------
% Parameters
% ----------------------------------------------------------
a1 = 1.2272;   % [cm2] Area of outlet pipe (Tank 1)
a2 = 1.2272;   % [cm2] Area of outlet pipe (Tank 2)
a3 = 1.2272;   % [cm2] Area of outlet pipe (Tank 3)
a4 = 1.2272;   % [cm2] Area of outlet pipe (Tank 4)
A1 = 380.1327; % [cm2] Cross sectional area (Tank 1)
A2 = 380.1327; % [cm2] Cross sectional area (Tank 2)
A3 = 380.1327; % [cm2] Cross sectional area (Tank 3)
A4 = 380.1327; % [cm2] Cross sectional area (Tank 4)
g = 981;       % [cm/s2] Gravity
gamma1 = 0.6;  % Flow distribution constant (valve 1)
gamma2 = 0.7;  % Flow distribution constant (valve 2)
rho = 1.0;     % [g/cm^3] Density of water

% Parameter-vector
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho]; 
at = p(1:4);   % [cm2] Area of outlet pipe
At = p(5:8);   % [cm2] Cross sectional area

% -----------------------------------------------------------
% Simulation scenario
% -----------------------------------------------------------
t0 = 0.0;           % [s] Initial time
tf= 20*60;          % [s] End time
N = 200;             % Number of steps 
dt = (tf-t0)/N;     % [s] interval between each step
t = t0:dt:tf;       % [s] time-vector
m10 = 0;            % [g] Liquid mass in tank 1 at time t0
m20 = 0;            % [g] Liquid mass in tank 2 at time t0
m30 = 0;            % [g] Liquid mass in tank 3 at time t0
m40 = 0;            % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;         % [cm3/s] Flow rate from pump 1
F2_0 = 300;         % [cm3/s] Flow rate from pump 2
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [0; 0;];                 % [cm3/s] Disturbance variables at t0
F_0 = [u0',d0'];              % [cm3/s] MV & DV at t0
u = u0.*ones(2, length(t));
d = d0.*ones(2, length(t));

% Solve ODE for this step
[T, X] = ode15s(@FourTankProcess, [t0:dt:tf], x0, [], u0, d0, p);

[y] = sensor_wo_noise(X, At, rho);
[z] = output(X, At, rho);
plots(t,u,y)

%%
% -------------------------- 2.2 -----------------------------------------%

d1 = 2*randn(1,length(t));             
d2 = 2*randn(1,length(t)); 
d_norm = [d1; d2];
u = u0.*ones(2, length(t));

[T, X, D_norm, U, x_norm] = discrete_fourtankProcess(x0, t, u, d_norm, p);

R = [1^2 0 0 0; 0 1^2 0 0; 0 0 0.5^2 0; 0 0 0 0.5^2];     % Covariance for disturbances in F3 and F4

[y_norm] = sensor_plus_noise(x_norm, At, rho, R);
[z_norm] = output(x_norm, At, rho);
plots(t,u,y_norm)

%---------Figure showing disturbance d1 and d2 in dicrete time-----------
figure(5)
plot(t/60, d_norm, 'LineWidth', 1);
xlabel('\textbf{t [min]}', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('[cm^{3}/s]', 'FontSize', 10);
xlim([0 t(end)/60]);
legend('F_3', 'F_4', 'Location', 'best');
title(['Disturbance following normal distribution (Discrete time with steps of ', ...
    num2str(dt), ' s)'], 'FontSize', 10);
%-------------------------------------------------------------------------

%---------Figure showing disturbance d1 and d2 in Continuous-time--------
figure(6)
plot(T/60, D_norm, 'LineWidth', 1);
xlabel('\textbf{t [min]}', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('[cm^{3}/s]', 'FontSize', 10);
xlim([0 T(end)/60]);
legend('F_3', 'F_4', 'Location', 'best');
title('Disturbance following normal distribution (Continuous-time)', ...
    'FontSize', 10);

%%
% -------------------------- 2.3 -----------------------------------------%

Ns = length(d0); % Number of realizations
seed = 100;
[W,t,dW] = ScalarStdWienerProcess(tf,N,Ns,seed);

sigma = [2^2 0; 0 2^2];                             % Covariance for disturbances in F3 and F4
d_brownian = sigma*dW';  


[T, X, D_brownian, U, x_brownian] = discrete_fourtankProcess(x0, t, u, d_brownian, p);


[y_brownian] = sensor_plus_noise(x_brownian, At, rho, R);

[z_brownian] = output(x_brownian, At, rho);

plots(t,u,y_brownian)

%---------Figure showing disturbance d1 and d2 in dicrete time-----------
figure(9)
plot(t/60, d_brownian, 'LineWidth', 1);
xlabel('\textbf{t [min]}', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('[cm^{3}/s]', 'FontSize', 10);
xlim([0 t(end)/60]);
legend('F_3', 'F_4', 'Location', 'best');
title(['Disturbance following Brownian motion (Discrete time with steps of ', ...
    num2str(dt), ' s)'], 'FontSize', 10);
%-------------------------------------------------------------------------

%---------Figure showing disturbance d1 and d2 in Continuous-time--------
figure(10)
plot(T/60, D_brownian, 'LineWidth', 1);
xlabel('\textbf{t [min]}', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('[cm^{3}/s]', 'FontSize', 10);
xlim([0 T(end)/60]);
legend('F_3', 'F_4', 'Location', 'best');
title('Disturbance following Brownian motion (Continuous-time)', ...
    'FontSize', 10);

%%
% -------------------------- 2.4 -----------------------------------------%


d0 = [0; 0;];
d = d0.*ones(2, length(t));
u(1,20:40) = 400;
u(2,100:end) = 100;

[T, X, D, U, x_step] = discrete_fourtankProcess(x0, t, u, d, p);
[y_step1] = sensor_wo_noise(x_step, At, rho);
plots(t,u,y_step1)
[y_step2] = sensor_plus_noise(x_step, At, rho, R);
plots(t,u,y_step2)


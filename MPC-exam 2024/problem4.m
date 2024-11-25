clc, clear, close all
addpath("Functions");

% Problem 4
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
tf= 160*60;                  % [s] End time
dt = 10;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
m10 = 0;                    % [g] Liquid mass in tank 1 at time t0
m20 = 0;                    % [g] Liquid mass in tank 2 at time t0
m30 = 0;                    % [g] Liquid mass in tank 3 at time t0
m40 = 0;                    % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
x0 = [m10; m20; m30; m40];  % [g] Start values 
u0 = [F1_0; F2_0];          % [cm3/s] Manipulated variables 
d0 = [0; 0;];               % [cm3/s] Disturbance variables at t0
d = d0.*ones(2, length(t));
u = u0.*ones(2, length(t));


%%  -------------------- 4.1 ---------------------------------------------
u_implement1 = find(t==40*60);
u(1,u_implement1:end) = u0(1)*(1-0.1);

u_implement2 = find(t==80*60);
u(1,u_implement2:end) = u(1,end)*(1+0.25);

u_implement3 = find(t==120*60);
u(1,u_implement3:end) = u(1,end)*(1-0.5);


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y] = sensor_wo_noise(x', At, rho);
[z] = output(x, At, rho);
plots(t,u,y')


u = u0.*ones(2, length(t));

u_implement1 = find(t==40*60);
u(2,u_implement1:end) = u0(2)*(1-0.1);

u_implement2 = find(t==80*60);
u(2,u_implement2:end) = u(2,end)*(1+0.25);

u_implement3 = find(t==120*60);
u(2,u_implement3:end) = u(2,end)*(1-0.5);


% Solve ODE for this step
[T2, X2, D2, U2, x2] = discrete_fourtankProcess(x0, t, u, d, p);

[y2] = sensor_wo_noise(x2', At, rho);
[z2] = output(x2, At, rho);
plots(t,u,y2')




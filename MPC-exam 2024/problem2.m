clc, clear

addpath("Functions");

%% Problem 2
% ----------------------------------------------------------
% Parameters
% ----------------------------------------------------------
a1 = 1.2272;    %[cm2] Area of outlet pipe (Tank 1)
a2 = 1.2272;    %[cm2] Area of outlet pipe (Tank 2)
a3 = 1.2272;    %[cm2] Area of outlet pipe (Tank 3)
a4 = 1.2272;    %[cm2] Area of outlet pipe (Tank 4)
A1 = 380.1327;  %[cm2] Cross sectional area (Tank 1)
A2 = 380.1327;  %[cm2] Cross sectional area (Tank 2)
A3 = 380.1327;  %[cm2] Cross sectional area (Tank 3)
A4 = 380.1327;  %[cm2] Cross sectional area (Tank 4)
g = 981;        %[cm/s2] Gravity
gamma1 = 0.6;  % Flow distribution constant (valve 1)
gamma2 = 0.7;  % Flow distribution constant (valve 2)
rho = 1.0;      % [g/cm^3] Density of water

p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho]; % Parameters
at = p(1:4);
At = p(5:8);

% -----------------------------------------------------------
% Simulation scenario
% -----------------------------------------------------------
t0 =    0.0;               % [s] Initial time
% t_int = 10                 % [s] interval time
tf= 10*60;                 % [s] End time
m10 = 0;                   % [g] Liquid mass in tank 1 at time t0
m20 = 0;                   % [g] Liquid mass in tank 2 at time t0
m30 = 0;                   % [g] Liquid mass in tank 3 at time t0
m40 = 0;                   % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                % [cm3/s] Flow rate from pump 1
F2_0 = 300;                % [cm3/s] Flow rate from pump 2

 
x0 = [m10; m20; m30; m40];    % Start values 
u0 = [F1_0; F2_0];            % Manipulated variables 
d0 = [0; 0;];                 % [cm3/s] Disturbance variables introduced at time 0
F_0 = [u0',d0'];              % [cm3/s] Manipulated variables and Disturbance variables introduced at time 0

% -------------------------- 2.1 -----------------------------------------%
% Solve ODE for this step
[t, X] = ode15s(@FourTankProcess, [t0 tf], x0, [], u0, d0, p);

Y = mass_to_height(X, [At; rho]);  % Liquid heights in all tanks

Z = mass_to_height(X(:,1:2), [At(1:2); rho]);  % Liquid heights in all tanks

% -------------------------- 2.2 -----------------------------------------%

d1=20*randn(length(t),1);             
d2=20*randn(length(t),1);              
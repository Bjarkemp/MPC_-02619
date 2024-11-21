clc; clear all; close all;
set(0,'DefaultTextInterpreter','latex'); set(0,'DefaultAxesFontSize',20); set(0,'DefaultLineLineWidth', 2);
objectType = 'Stair';
propertyName = 'LineWidth';
set(groot, ['Default', objectType, propertyName], 2)

%% Linearize 
syms m1 m2 m3 m4 F1 F2 F3 F4 h1 h2 h3 h4 gamma1 gamma2 A1 A2 A3 A4 g rho a1 a2 a3 a4

eq1 = rho*((gamma1*F1) + (a3.*sqrt(2*g*m3/(rho*A3))) - (a1.*sqrt(2*g*m1/(rho*A1)))); % Mass balance Tank 1
eq2 = rho*((gamma2*F2) + (a4.*sqrt(2*g*m4/(rho*A4))) - (a2.*sqrt(2*g*m2/(rho*A2)))); % Mass balance Tank 2
eq3 = rho*(((1-gamma2)*F2) - (a3.*sqrt(2*g*m3/(rho*A3))) + F3); % Mass balance Tank 3
eq4 = rho*(((1-gamma1)*F1) - (a4.*sqrt(2*g*m4/(rho*A4))) + F4); % Mass balance Tank 4
    
eq1_out = m1/(rho*A1);
eq2_out = m2/(rho*A2);
eq1_out_h = h1;
eq2_out_h = h2;

% Linearization
A_alg = jacobian([eq1 eq2 eq3 eq4] , [m1 m2 m3 m4])
B_alg = jacobian([eq1 eq2 eq3 eq4] , [F1 F2])
E_alg = jacobian([eq1 eq2 eq3 eq4] , [F3 F4])
C_alg = jacobian([eq1_out eq2_out] , [m1 m2 m3 m4])
C_h_alg = jacobian([eq1_out_h eq2_out_h] , [h1 h2 h3 h4])
D_alg = jacobian([eq1_out eq2_out] , [F1 F2])


%% Parameters
% --------------------------------------------------------------
a1 = 1.2272; %[cm2] Area of outlet pipe 1
a2 = 1.2272; %[cm2] Area of outlet pipe 2
a3 = 1.2272; %[cm2] Area of outlet pipe 3
a4 = 1.2272; %[cm2] Area of outlet pipe 4
A1 = 380.1327; %[cm2] Cross sectional area of tank 1
A2 = 380.1327; %[cm2] Cross sectional area of tank 2
A3 = 380.1327; %[cm2] Cross sectional area of tank 3
A4 = 380.1327; %[cm2] Cross sectional area of tank 4
gamma1 = 0.58; % Flow distribution constant. Valve 1
gamma2 = 0.68; % Flow distribution constant. Valve 2
g = 981; %[cm/s2] The acceleration of gravity
rho = 1.00; %[g/cm3] Density of water
p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];
% --------------------------------------------------------------


% Initial values
% --------------------------------------------------------------
t0 = 0; % [s] Initial time
t_f = 20*60; % [s] Final time
m10 = 0.0; % [g] Liquid mass in tank 1 at time t0
m20 = 0.0; % [g] Liquid mass in tank 2 at time t0
m30 = 0.0; % [g] Liquid mass in tank 3 at time t0
m40 = 0.0; % [g] Liquid mass in tank 4 at time t0
x0 = [m10; m20; m30; m40];

%% Steady state value
u = [300 ; 300];
d = [250 ; 250];
[T,X] = ode15s(@ModifiedFourTankSystem,[t0 t_f],x0,[],u,d,p); % REMEBMER ODE OPTIONS SHOULD []!!!

% Steady state value
x_s = X(end,:); % [g]

%% Linear model numerical
% Subs steady state values
m1 = x_s(1); m2 = x_s(2); m3 = x_s(3); m4 = x_s(4);
A_c = double(subs(A_alg));
B_c = double(subs(B_alg));
E_c = double(subs(E_alg));
C_c = double(subs(C_alg));
C_z = double(subs(C_h_alg));
D_c = double(subs(D_alg));

%% Discretize
Ts = 15; % Given in assignment
mat = expm([A_c , B_c  ,E_c ; zeros(size(B_c,2),size(A_c,1)) , zeros(size(B_c,2)) , zeros(size(E_c,2)) ; zeros(size(E_c,2),size(A_c,1)) , zeros(size(B_c,2)) , zeros(size(E_c,2))].*Ts);
A_k = mat(1:4,1:4); B_k = mat(1:4,5:6); E_k = mat(1:4,7:8);

save('Parameters')

%% Function
function xdot = ModifiedFourTankSystem(t,x0,u,d,p)
    % MODIFIEDFOURTANKSYSTEM Model dx/dt = f(t,x,u,d,p) for Modified 4-tank System
    %
    % This function implements a differential equation model for the
    % modified 4-tank system.
    %
    % Syntax: xdot = ModifiedFourTankSystem(t,x,u,d,p)
    % Unpack states, MVs, and parameters
    m = x0; % Mass of liquid in each tank [g]
    F = [u; d]; % Flow rates [cm3/s]
    a = p(1:4,1); % Pipe cross sectional areas [cm2]
    A = p(5:8,1); % Tank cross sectional areas [cm2]
    gamma = p(9:10,1); % Valve positions [-]
    g = p(11,1); % Acceleration of gravity [cm/s2]
    rho = p(12,1); % Density of water [g/cm3]
    % Inflows
    qin = zeros(4,1);
    qin(1,1) = gamma(1)*F(1); % Inflow from valve 1 to tank 1 [cm3/s]
    qin(2,1) = gamma(2)*F(2); % Inflow from valve 2 to tank 2 [cm3/s]
    qin(3,1) = (1-gamma(2))*F(2); % Inflow from valve 2 to tank 3 [cm3/s]
    qin(4,1) = (1-gamma(1))*F(1); % Inflow from valve 1 to tank 4 [cm3/s]
    % Outflows
    h = m./(rho*A); % Liquid level in each tank [cm]
    qout = a.*sqrt(2*g*h); % Outflow from each tank [cm3/s]
    % Differential equations
    xdot = zeros(4,1);
    xdot(1,1) = rho*(qin(1,1)+qout(3,1)-qout(1,1)); % Mass balance Tank 1
    xdot(2,1) = rho*(qin(2,1)+qout(4,1)-qout(2,1)); % Mass balance Tank 2
    xdot(3,1) = rho*(qin(3,1)-qout(3,1)+F(3)); % Mass balance Tank 3
    xdot(4,1) = rho*(qin(4,1)-qout(4,1)+F(4)); % Mass balance Tank 4
end


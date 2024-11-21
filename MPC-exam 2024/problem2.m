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
tf= 10*60;          % [s] End time
N = 50;             % Number of steps 
t_int = (tf-t0)/N;  % [s] interval between each step
t = t0:t_int:tf;    % [s] time-vector
m10 = 0;            % [g] Liquid mass in tank 1 at time t0
m20 = 0;            % [g] Liquid mass in tank 2 at time t0
m30 = 0;            % [g] Liquid mass in tank 3 at time t0
m40 = 0;            % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;         % [cm3/s] Flow rate from pump 1
F2_0 = 300;         % [cm3/s] Flow rate from pump 2

 
x0 = [m10; m20; m30; m40];    % Start values 
u0 = [F1_0; F2_0];            % Manipulated variables 
d0 = [0; 0;];                 % [cm3/s] Disturbance variables at t0
F_0 = [u0',d0'];              % [cm3/s] MV & DV at t0

% -------------------------- 2.1 -----------------------------------------%
% Solve ODE for this step
[T, X] = ode15s(@FourTankProcess, [t0 tf], x0, [], u0, d0, p);

Y = mass_to_height(X, [At; rho]);              % Heights in all tanks
Z = mass_to_height(X(:,1:2), [At(1:2); rho]);  % Heights in all tanks
u = u0.*ones(2, length(T));
plots(T,u,Y)

%%
% -------------------------- 2.2 -----------------------------------------%

d1=20*randn(1,length(t));             
d2=20*randn(1,length(t)); 
d=[d1; d2];
u = u0.*ones(2, length(t));

% Initialising compiled masses and time
x_sample = x0';     % [g] Mass at each sample 
X =[];              % [g] Mass at all times
T = [];             % [s] Same time period, but with much more time steps which does not represent time between samples
D = [];
U = [];
seed = 1000;
rng(seed)           % (du kan vælge et andet heltal som seed)

% Discritized solution of four tank system.
for k=1:length(t)-1
    [Tk, Xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], x0, [], u(:,k), d(:,k), p);

    x0=Xk(end,:);                                  % The last value of Xk is start value for next iteration
    x_sample = [x_sample; Xk(end,:)];              % Compiling mass for each sample after every iteration

    T = [T; Tk];                                   % Compiling time after every iteration
    X = [X; Xk];                                   % Compiling  mass after every iteration
    D = [D; (d(:,k))'.*ones(length(Xk),2)];
    U = [U; (u(:,k))'.*ones(length(Xk),2)]; 
end

y_sample = mass_to_height(x_sample, [At; rho]);              % Liquid heights in all tanks
plots(t,u,y_sample)

figure(5)
plot(t/60, d, 'LineWidth', 1);
xlabel('\textbf{t [min]}', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('[cm^{3}/s]', 'FontSize', 10);
xlim([0 t(end)/60]);
legend('d_1', 'd_2', 'Location', 'best');
title(['Disturbance in F_3 and F_4 (Discrete time with steps of ', num2str(t_int), ' s)'], 'FontSize', 10);

figure(6)
plot(T/60, D, 'LineWidth', 1);
xlabel('\textbf{t [min]}', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('[cm^{3}/s]', 'FontSize', 10);
xlim([0 T(end)/60]);
legend('d_1', 'd_2', 'Location', 'best');
title('Disturbance in F_3 and F_4 to tank 3 and 4 (Continuous-time)', 'FontSize', 10);
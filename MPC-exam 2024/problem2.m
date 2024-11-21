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

% -------------------------- 2.1 -----------------------------------------%
% Solve ODE for this step
[T, X] = ode15s(@FourTankProcess, [t0 tf], x0, [], u0, d0, p);

Y = mass_to_height(X, [At; rho]);              % Heights in all tanks
Z = mass_to_height(X(:,1:2), [At(1:2); rho]);  % Heights in all tanks
u = u0.*ones(2, length(T));
plots(T,u,Y)

%%
% -------------------------- 2.2 -----------------------------------------%

d1=5*randn(1,length(t));             
d2=5*randn(1,length(t)); 
d=[d1; d2];
u = u0.*ones(2, length(t));

% Initialising compiled masses and time
x_discrete = x0';     % [g] Mass at each dicrete time 
T = [];             % [s] Continuous-time
X =[];              % [g] Mass at all at Continuous-time
D = [];             % [cm3/s] Disturbance for Continuous-time
U = [];             % [cm3/s] Manipulated variables

seed = 1000;
rng(seed)           % Seed for randon number-generation


% Discritized solution of four tank system.
for k=1:length(t)-1
    [Tk, Xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], x0, ...
        [], u(:,k), d(:,k), p);

    x0=Xk(end,:);                     % Xk(end:) is Discrete time sample
    x_discrete = [x_discrete; Xk(end,:)]; % [g] Mass for Discrete time

    T = [T; Tk];                      
    X = [X; Xk];                      
    D = [D; (d(:,k))'.*ones(length(Xk),2)];
    U = [U; (u(:,k))'.*ones(length(Xk),2)]; 
end


R = [1^2 0 0 0; 0 1^2 0 0; 0 0 0.5^2 0; 0 0 0 0.5^2];     % Covariance for disturbances in F3 and F4
Lr = chol(R,'lower');                                     % Cholesky-dekomposition. It just gives me the standard deviation instead of variance.
v = (Lr*randn(length(x0),length(t)))';                    % Measurement noise. Follows normal distribution with mean=0 and has st.dev of Lr
v(1,:)=zeros(1,length(x0));                               % no measurement noise at t=0

y_sample = mass_to_height(x_discrete, [At; rho])+v; % Heights in all tanks

z = mass_to_height(x_discrete(:,1:2), [At(1:2); rho]);  % Heights in all tanks

plots(t,u,y_sample)

%---------Figure showing disturbance d1 and d2 in dicrete time-----------
figure(5)
plot(t/60, d, 'LineWidth', 1);
xlabel('\textbf{t [min]}', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('[cm^{3}/s]', 'FontSize', 10);
xlim([0 t(end)/60]);
legend('d_1', 'd_2', 'Location', 'best');
title(['Disturbance in F_3 and F_4 (Discrete time with steps of ', ...
    num2str(dt), ' s)'], 'FontSize', 10);
%-------------------------------------------------------------------------

%---------Figure showing disturbance d1 and d2 in Continuous-time--------
figure(6)
plot(T/60, D, 'LineWidth', 1);
xlabel('\textbf{t [min]}', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('[cm^{3}/s]', 'FontSize', 10);
xlim([0 T(end)/60]);
legend('d_1', 'd_2', 'Location', 'best');
title('Disturbance in F_3 and F_4 to tank 3 and 4 (Continuous-time)', ...
    'FontSize', 10);
%-------------------------------------------------------------------------
%%
% -------------------------- 2.3 -----------------------------------------%

Ns = length(x0); % Number of realizations
T = tf; % Final time
N = length(t)-1; % Number of time steps
seed = 1*randn+20; % Seed for reproducibility
rng(seed); % seed for my random generation
dt = T/N; % difference in time for each interval
[W,tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed);
[Wmean,sW,Wmeanp2s,Wmeanm2s]=ScalarSampleMeanStdVar(W);
clc, clear, %close all
addpath("Functions");

% Problem 3
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
tf= 3*60;                  % [s] End time
dt = 1;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
m10 = 17612.0123865384;                    % [g] Liquid mass in tank 1 at time t0
m20 = 29640.6694949484;                    % [g] Liquid mass in tank 2 at time t0
m30 = 4644.21948249808;                    % [g] Liquid mass in tank 3 at time t0
m40 = 9378.49308238605; 
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 =100;
F4_0 =150;
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];           % [cm3/s] Disturbance variables at t0
y0 = mass_to_height(x0,At,rho);
d = d0.*ones(2, length(t));
u = u0.*ones(2, length(t));


xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);    % Løser differentialligningssystemet
% (QuadrupleTankProcess) og solver hvad x er når hældningen er 0. Dvs. at
% den beregner hvad masserne er når der opnås steady state i tankene.


%%


% Modelling the disturbance as Brownian motion
Ns = length(d0); % Number of realizations
seed = 100;
[W,t,dW] = ScalarStdWienerProcess(tf,N,Ns,seed);
sigma = [0.0^2 0; 0 0.0^2];                             % Covariance for disturbances in F3 and F4
d = d0 + sigma*dW'; 


% Step changes in manipulated variables
% u(2,30:end) = u(2,30)*1.5;
% u(2,300:end) = u(2,300)*0.8;
% u(1,900:end) = u(1,900)*1.15;

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;           % Covariance for process noise

[T, X, D, U, x_sample] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

% Deviation variables
udev = u - u0;
xdev = x_sample' - x0;
ddev = d - d0;


% --------------------------------------------------------------
% Linearization
% --------------------------------------------------------------
ys = mass_to_height(xs,At,rho);
hs = ys;
T = (At./at).*sqrt(2*hs/g);
A=[-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B=[rho*gamma1 0;0 rho*gamma2; 0 rho*(1-gamma2); rho*(1-gamma1) 0];
C=diag(1./(rho*At));
Cz=C(1:2,:);
G = [0 0; 0 0; rho 0; 0 rho];

% ZOH Discretization of Linear System
[Ad, Bd, Gd] = c2dzoh(A, B, G, dt);


% Kalman filter
[x_hat, x_phat] = kalman_filter(t, xdev, udev, ddev, At, rho, R, Q, Ad, Bd, Gd, C);

% From deviation variables to 
x = xdev + x0;
x_hat = x_hat + x0;

figure(1)
for i = 1:4
    subplot(2,2,i)
    plot(t/60, x(i,:),t/60, x_hat(i,:), 'LineWidth', 2);
    grid on;
    xlabel('t [min]', 'FontSize', 12);
    ylabel('m [g]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    legend('Mass in tank', 'Kalman filter', 'Location', 'best');
    title(['Mass in tank ', num2str(i)], 'FontSize', 10);
end



%%


A_aug = [A, G; zeros(size(G, 2), size(A, 1)), eye(size(G, 2))];
B_aug = [B; zeros(size(G, 2), size(B, 2))];
G_aug = [G; zeros(size(G, 2), size(G, 2))];
C_aug = [C, zeros(size(C, 1), size(G, 2))];
Q_aug = Q;  % Diagonal kovariansmatrix for forstyrrelser
G_aug = [G zeros(4,2); zeros(2,2) eye(2,2)];



% ZOH Discretization of Linear System
[Ad_aug, Bd_aug, Gd_aug] = c2dzoh(A_aug, B_aug, G_aug, dt);
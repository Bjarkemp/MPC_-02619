clc, clear, close all
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
tf= 50*60;                  % [s] End time
dt = 10;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
m10 = 17612.0123865358;     % [g] Liquid mass in tank 1 at time t0
m20 = 29640.6694949367;     % [g] Liquid mass in tank 2 at time t0
m30 = 4644.21948249842;     % [g] Liquid mass in tank 3 at time t0
m40 = 9378.49308238594;     % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 =100;
F4_0 =150;
x0 = [m10; m20; m30; m40];  % [g] Start values 
u0 = [F1_0; F2_0];          % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];         % [cm3/s] Disturbance variables at t0
y0 = mass_to_height(x0,At,rho);
d = d0.*ones(2, length(t));
u = u0.*ones(2, length(t));


xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);    % Løser differentialligningssystemet
% (QuadrupleTankProcess) og solver hvad x er når hældningen er 0. Dvs. at
% den beregner hvad masserne er når der opnås steady state i tankene.


%% Linear Kalman filter (Not Augmented)

% Modelling the disturbance as Brownian motion
Ns = length(d0); % Number of realizations
seed = 10;
[W,t,dW] = ScalarStdWienerProcess(tf,N,Ns,seed);
sigma = [2^2 0; 0 2^2];                             % Covariance for disturbances in F3 and F4
% d = d0 + sigma*dW'; 
% d = [130*ones(1,length(t));190*ones(1,length(t))];

% Step changes in manipulated variables
u(2,1:end) = u(2,1)*0.5;
u(1,1:end) = u(1,1)*1.5;
% d(2,250:end) = d(2,250:end)+100;


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
Tl = (At./at).*sqrt(2*hs/g);
A=[-1/Tl(1) 0 1/Tl(3) 0;0 -1/Tl(2) 0 1/Tl(4);0 0 -1/Tl(3) 0;0 0 0 -1/Tl(4)];
B=[rho*gamma1 0;0 rho*gamma2; 0 rho*(1-gamma2); rho*(1-gamma1) 0];
C=diag(1./(rho*At));
Cz=C(1:2,:);
E = [0 0; 0 0; rho 0; 0 rho];
Gw = eye(4);

%-------------------------------------------------------------------------
% Construct the matrix for exponential calculation
M = expm([-A E*E' ; zeros(size(A)) A']*dt);

% Extract submatrices
phi_12 = M(1:4,5:8);
phi_22 = M(5:8,5:8);

% Compute process noise covariance matrix Q
Q2 = phi_22' * phi_12;
%-------------------------------------------------------------------------
% ZOH Discretization of Linear System
[Ad, Bd, Ed] = c2dzoh(A, B, E, dt);


% Kalman filter
[x_hat1_dyn, x_phat1_dyn] = kalman_filter_dynamic(t, xdev, udev, ddev, At, rho, R, Q, Ad, Bd, Ed, C);
[x_hat1_sta, x_phat1_sta] = kalman_filter_static(t, xdev, udev, ddev, At, rho, R, Q, Ad, Bd, Ed, C);

% From deviation variables to 
x = xdev + x0;
xhat1_dyn = x_hat1_dyn + x0;
xhat1_sta = x_hat1_sta + x0;

y = sensor_plus_noise(x',At,rho,R);
yhat1_dyn = mass_to_height(xhat1_dyn(1:4,:),At,rho);
yhat1_sta = mass_to_height(xhat1_sta(1:4,:),At,rho);


figure(1)
for i = 1:4
    subplot(2,2,i)
    plot(t/60, y(i,:),'b', 'LineWidth', 2); 
    hold on;
    plot(t/60, yhat1_dyn(i,:),'r', 'LineWidth', 1);
    hold on;
    plot(t/60, yhat1_sta(i,:),'g', 'LineWidth', 1);
    hold off;
    grid on;
    ylabel('height [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    legend('Measured height', 'Dynamic Kalman filter', 'Static Kalman filter', 'Location', 'best');
    title(['Tank ', num2str(i)], 'FontSize', 10);
end
sgtitle('Not-Augmented Linear Kalman filter', 'FontSize', 14, 'FontWeight', 'bold');

%% Augmented Linear Kalman filter 

% Diskretisering af det kontinuerte system
[Ad, Bd, Ed] = c2dzoh(A, B, E, dt);
% ZOH Discretization of Linear System
[Ad, Bd, Gwd] = c2dzoh(A, B, Gw, dt);

% Augmentering efter diskretisering
Ad_aug = [Ad, Ed; zeros(size(Ed, 2), size(Ad, 1)), eye(size(Ed, 2))];
Bd_aug = [Bd; zeros(size(Ed, 2), size(Bd, 2))];
Ed_aug = [Ed; zeros(size(Ed, 2), size(Ed, 2))];
C_aug = [C, zeros(size(C, 1), size(Ed, 2))];
Q_aug = [Q, zeros(size(Q, 1), size(Ed, 2)); zeros(size(Ed, 2), size(Q, 2)), eye(size(Ed, 2))];
Gw_aug =  [Gw, zeros(4,2); zeros(2,4) eye(2,2)];

[x_hat2_dyn, x_phat2_dyn] = kalman_filter_aug_dynamic(t, xdev, udev, ddev, At, rho, R, Q_aug, Ad_aug, Bd_aug, Ed_aug, Gw_aug, C_aug);
[x_hat2_sta, x_phat2_sta] = kalman_filter_aug_static(t, xdev, udev, ddev, At, rho, R, Q_aug, Ad_aug, Bd_aug, Ed_aug, Gw_aug, C_aug);

[x_hat2_dyn_pre, x_phat2_dyn_pre] = kalman_filter_aug_dynamic_pred(t(1), xdev(:,1), udev(:,1), At, rho, R, Q_aug, Ad_aug, Bd_aug, Gw_aug, C_aug,300);

% From deviation variables to 
xhat2_dyn = x_hat2_dyn(1:4,:) + x0;
xhat2_sta = x_hat2_sta(1:4,:) + x0;
xhat2_dyn_pre = x_phat2_dyn_pre(1:4,:) + x0;

yhat2_dyn = mass_to_height(xhat2_dyn(1:4,:),At,rho);
yhat2_sta = mass_to_height(xhat2_sta(1:4,:),At,rho);
yhat2_dyn_pre = mass_to_height(xhat2_dyn_pre(1:4,:),At,rho);

figure(2)
for i = 1:4
    subplot(2,2,i)
    plot(t/60, y(i,:),'b', 'LineWidth', 2); 
    hold on;
    plot(t/60, yhat2_dyn(i,:),'r', 'LineWidth', 1);
    hold on;
    plot(t/60, yhat2_sta(i,:),'g', 'LineWidth', 1);
    hold on;
    plot(t(1:300)/60, yhat2_dyn_pre(i,:),'k', 'LineWidth', 1);
    hold off;
    grid on;
    xlabel('t [min]', 'FontSize', 12);
    ylabel('height [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    legend('Measured height', 'Dynamic Kalman filter', 'Static Kalman filter', 'Predicting Dynamic Kalman filter', 'Location', 'best');
    title(['Tank ', num2str(i)], 'FontSize', 10);
end
sgtitle('Augmented Linear Kalman filter', 'FontSize', 14, 'FontWeight', 'bold');
%% Augmented non-Linear Kalman filter 


[x_hat3_dyn, x_phat3_dyn] = nonlinear_kalman_filter_aug_dynamic(t, xdev, udev, ddev, x0, u, d, At, rho, R, Q_aug, Ad_aug, Gw_aug, C_aug, p);
[x_hat3_sta, x_phat3_sta] = nonlinear_kalman_filter_aug_static(t, xdev, udev, ddev,  x0, u, d, At, rho, R, Q_aug, Ad_aug, Gw_aug, C_aug, p);


% From deviation variables to 
xhat3_dyn = x_hat3_dyn(1:4,:) + x0;
xhat3_sta = x_hat3_sta(1:4,:) + x0;

yhat3_dyn = mass_to_height(xhat3_dyn(1:4,:),At,rho);
yhat3_sta = mass_to_height(xhat3_sta(1:4,:),At,rho);

figure(3)
for i = 1:4
    subplot(2,2,i)
    plot(t/60, y(i,:),'b', 'LineWidth', 2); 
    hold on;
    plot(t/60, yhat3_dyn(i,:),'r', 'LineWidth', 1);
    hold on;
    plot(t/60, yhat3_sta(i,:),'g', 'LineWidth', 1);
    hold off;
    grid on;
    xlabel('t [min]', 'FontSize', 12);
    ylabel('height [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    legend('Measured height', 'Dynamic Kalman filter', 'Static Kalman filter', 'Location', 'best');
    title(['Tank ', num2str(i)], 'FontSize', 10);
end
sgtitle('Augmented Non-linear Kalman filter', 'FontSize', 14, 'FontWeight', 'bold');







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
tf= 50*60;                  % [s] End time
dt = 10;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
m10 = 1000;     % [g] Liquid mass in tank 1 at time t0
m20 = 2000;     % [g] Liquid mass in tank 2 at time t0
m30 = 3000;     % [g] Liquid mass in tank 3 at time t0
m40 = 8000;     % [g] Liquid mass in tank 4 at time t0
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


%%

% Modelling the disturbance as Brownian motion
Ns = length(d0); % Number of realizations
seed = 10;
[W,t,dW] = ScalarStdWienerProcess(tf,N,Ns,seed);
sigma = [2^2 0; 0 2^2];                             % Covariance for disturbances in F3 and F4
d = d0 + sigma*dW'; 


% Step changes in manipulated variables
u(2,50:end) = u(2,59)*0.5;
u(1,200:end) = u(1,200)*1.25;
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
Tl = (At./at).*sqrt(2*hs/g);
A=[-1/Tl(1) 0 1/Tl(3) 0;0 -1/Tl(2) 0 1/Tl(4);0 0 -1/Tl(3) 0;0 0 0 -1/Tl(4)];
B=[rho*gamma1 0;0 rho*gamma2; 0 rho*(1-gamma2); rho*(1-gamma1) 0];
C=diag(1./(rho*At));
Cz=C(1:2,:);
G = [0 0; 0 0; rho 0; 0 rho];
Gw = eye(4);

% ZOH Discretization of Linear System
[Ad, Bd, Gd] = c2dzoh(A, B, G, dt);


% Kalman filter
[x_hat1, x_phat1] = kalman_filter_dynamic(t, xdev, udev, ddev, At, rho, R, Q, Ad, Bd, Gd, C);

% From deviation variables to 
x = xdev + x0;
xhat1 = x_hat1 + x0;

y1 = sensor_plus_noise(x',At,rho,R);
yhat1 = mass_to_height(xhat1(1:4,:),At,rho);
z1 = output(xhat1(1:4,:)',At,rho);

figure(1)
for i = 1:4
    subplot(2,2,i)
    plot(t/60, y1(i,:),'b', 'LineWidth', 2); 
    hold on;
    plot(t/60, yhat1(i,:),'r', 'LineWidth', 1);
    hold off;
    grid on;
    xlabel('t [min]', 'FontSize', 12);
    ylabel('h [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    legend('Height in tank', 'Kalman filter', 'Location', 'best');
    title(['Height in tank ', num2str(i)], 'FontSize', 10);
end
sgtitle('Not-Augmented Linear Kalman filter', 'FontSize', 14, 'FontWeight', 'bold');

%%

% Diskretisering af det kontinuerte system
[Ad, Bd, Gd] = c2dzoh(A, B, G, dt);
% ZOH Discretization of Linear System
[Ad, Bd, Gwd] = c2dzoh(A, B, Gw, dt);

% Augmentering efter diskretisering
Ad_aug = [Ad, Gd; zeros(size(Gd, 2), size(Ad, 1)), eye(size(Gd, 2))];
Bd_aug = [Bd; zeros(size(Gd, 2), size(Bd, 2))];
Gd_aug = [Gd; zeros(size(Gd, 2), size(Gd, 2))];
C_aug = [C, zeros(size(C, 1), size(Gd, 2))];
Q_aug = [Q, zeros(size(Q, 1), size(Gd, 2)); zeros(size(Gd, 2), size(Q, 2)), eye(size(Gd, 2))];
Gw_aug =  [Gwd, zeros(4,2); zeros(2,6)];

[x_hat2, x_phat2] = kalman_filter_aug_dynamic(t, xdev, udev, ddev, At, rho, R, Q_aug, Ad_aug, Bd_aug, Gd_aug, Gw_aug, C_aug);


% From deviation variables to 
x = xdev + x0;
xhat2 = x_hat2(1:4,:) + x0;

y2 = sensor_plus_noise(x',At,rho,R);
yhat2 = mass_to_height(xhat2(1:4,:),At,rho);
z2 = output(xhat2(1:4,:)',At,rho);

figure(2)
for i = 1:4
    subplot(2,2,i)
    plot(t/60, y2(i,:),'b', 'LineWidth', 2); 
    hold on;
    plot(t/60, yhat2(i,:),'r', 'LineWidth', 1);
    hold off;
    grid on;
    xlabel('t [min]', 'FontSize', 12);
    ylabel('h [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    legend('Height in tank', 'Kalman filter', 'Location', 'best');
    title(['Height in tank ', num2str(i)], 'FontSize', 10);
end
sgtitle('Augmented Linear Kalman filter', 'FontSize', 14, 'FontWeight', 'bold');



[x_hat3, x_phat3] = nonlinear_kalman_filter_aug_dynamic(t, xdev, udev, ddev, x0, u, d, At, rho, R, Q_aug, Ad_aug, Gw_aug, C_aug, p);

% From deviation variables to 
x = xdev + x0;
xhat3 = x_hat3(1:4,:) + x0;

y3 = sensor_plus_noise(x',At,rho,R);
yhat3 = mass_to_height(xhat3(1:4,:),At,rho);
z3 = output(xhat3(1:4,:)',At,rho);

figure(3)
for i = 1:4
    subplot(2,2,i)
    plot(t(1:end)/60, y3(i,:),'b', 'LineWidth', 2); 
    hold on;
    plot(t(1:end)/60, yhat3(i,:),'r', 'LineWidth', 1);
    hold off;
    grid on;
    xlabel('t [min]', 'FontSize', 12);
    ylabel('h [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    legend('Height in tank', 'Kalman filter', 'Location', 'best');
    title(['Height in tank ', num2str(i)], 'FontSize', 10);
end
sgtitle('Augmented Non-linear Kalman filter', 'FontSize', 14, 'FontWeight', 'bold');

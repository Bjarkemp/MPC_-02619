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
tf= 60*10;                  % [s] End time
dt = 2;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
Ph = 5;                     % Prediction horizon
m10 = 0;                    % [g] Liquid mass in tank 1 at time t0
m20 = 0;                    % [g] Liquid mass in tank 2 at time t0
m30 = 0;                    % [g] Liquid mass in tank 3 at time t0
m40 = 0;                    % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 = 100;
F4_0 = 150;
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];           % [cm3/s] Disturbance variables at t0
d = d0.*ones(2, length(t));
[y0] = sensor_wo_noise(x0', At, rho);

%%  -------------------- 8.1 MPC function ------------------------------

%linearization
% Steady State
xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);
ys = sensor_wo_noise(xs,at,rho);
zs = sensor_wo_noise(xs,at,rho);

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;           % Covariance for process noise

[A, B, C, Cz, E, Gw] = system_matrices(xs, At, at, rho, gamma1, gamma2, g);

%ZOH Discretization of Linear System
%Stochastic Brownian
[Ad,Bd,Ed]=c2dzoh(A,B,E,dt);
D = zeros(2,4);
% sys = ss(Ad,[Bd,Gd],C(1:2,:),D);
sys = ss(Ad,Bd,C(1:2,:),D(:,1:2));


% Augmentation after discritization
Ad_aug = [Ad, Ed; zeros(size(Ed, 2), size(Ad, 1)), eye(size(Ed, 2))];
Bd_aug = [Bd; zeros(size(Ed, 2), size(Bd, 2))];
Ed_aug = [Ed; zeros(size(Ed, 2), size(Ed, 2))];
C_aug = [C, zeros(size(C, 1), size(Ed, 2))];
Q_aug = [Q, zeros(size(Q, 1), size(Ed, 2)); zeros(size(Ed, 2), size(Q, 2)), eye(size(Ed, 2))];
Gw_aug =  [Gw, zeros(4,2); zeros(2,4) eye(2,2)];

A = sys.A;
B = sys.B;
C = sys.C;

% Define MPC parameters
Qz = 10000 * eye(size(C, 1));  % Weight on output tracking
S = 1 * eye(size(B, 2));  % Weight on control effort

% Design MPC
MPC_sys = UnconstrainedMPCDesign(A, B, C, Qz, S, Ph);

h_sp = [20; 30];
hdev_sp = h_sp - y0(1:2);

% Rsp = 100*ones(2*Ph,1);
Rsp = [hdev_sp(1)*ones(1, Ph); hdev_sp(2)*ones(1, Ph)];  % To referencer for Tank 1 og Tank 2

% Ekstraher den relevante del af reference
ref_traj = reshape(Rsp(:, 1:min(1+Ph-1, end)), [], 1);


% Modelling the disturbance as Brownian motion
Ns = length(d0); % Number of realizations
seed = 10;
[W,t,dW] = ScalarStdWienerProcess(tf,N,Ns,seed);
sigma = [1^2 0; 0 1^2];                             % Covariance for disturbances in F3 and F4
d_k = d0 + sigma*dW';


inputs = {x0, u0, Rsp, d_k};

% Simulation
% [y, u] = MPC_Sim_Unconstrained(sys, MPC_sys, Q_aug, R, tf, dt, inputs, u0,t,At,rho, Ad_aug, Bd_aug, Ed_aug, Gw_aug, C_aug, Ph);


% Deviation variables
xdev = zeros(4,1);
udev(:,1) = zeros(2,1);
ddev = d_k - d(:,1);
%%

for i = 1:length(t)

% Kalman filter
%--------------------------------------------------------------------------
[x_hat, x_phat] = kalman_filter_aug_dynamic_pred(t(i), xdev(:,i), udev(:,i), ddev(:,i), At, rho, R, Q_aug, Ad_aug, Bd_aug, Ed_aug, Gw_aug, C_aug, Ph);
%--------------------------------------------------------------------------
xdev(:,i) = x_hat(1:4,1);
xdev(:,i+1) = x_phat(1:4,1);

x_mpc = [x_hat(1:4,1) x_phat(1:4,:)];

for j = 1:Ph
    % Beregn den lineære del af omkostningsfunktionen
    g = MPC_sys.M_x0 * x_mpc(:,j) + MPC_sys.M_r * ref_traj;
    
    % Løs QP-problemet for denne iteration
    u_current = qpsolver(MPC_sys.H, g, [], [], [], [], [], []);
    
    % Gem det første input fra løsningen
    u_mpc((j-1)*size(B, 2)+1:j*size(B, 2)) = u_current(1:size(B, 2));
end

    % Gem de opdaterede estimater og fremtidige prediktioner
    x_hat_log(:,i) = xdev(:,i) + x0;          % Gem sidste opdaterede estimat
    
    % Apply first control action
    % udev(:, i+1) = u_mpc(1:2); % First control input for the time step
    udev(:, i+1) = u_mpc(end-1:end);
    ulog = udev(:,1:end-1) + u0;
end

% Konverter til højder for visualisering
yhat_log = sensor_plus_noise(x_hat_log', At, rho, R);



figure(1)
for i = 1:2
    subplot(2,2,i)
    plot(t/60, yhat_log(i,:),'b', 'LineWidth', 2);
    hold on 
    grid on;
    ylabel('height [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    ylim([0 40]);
    % legend('Measured height', 'Dynamic Kalman filter', 'Static Kalman filter', 'Location', 'best');
    title(['Tank ', num2str(i)], 'FontSize', 10);
    hold on
end
for i = 1:2
    subplot(2,2,i+2)
    plot(t/60, ulog(i,:),'b', 'LineWidth', 2); 
    hold off;
    grid on;
    ylabel('[cm^3/s]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    % ylim([0 600]);
    % legend('Measured height', 'Dynamic Kalman filter', 'Static Kalman filter', 'Location', 'best');
    title(['F', num2str(i)], 'FontSize', 10);
    hold off
end
sgtitle('  ', 'FontSize', 14, 'FontWeight', 'bold');
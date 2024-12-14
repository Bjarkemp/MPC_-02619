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
tf= 60*20;                  % [s] End time
dt = 5;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
Ph = 15;                    % Prediction horizon
m10 = 17612.0123865384;                    % [g] Liquid mass in tank 1 at time t0
m20 = 29640.6694949484;                    % [g] Liquid mass in tank 2 at time t0
m30 = 4644.21948249808;                    % [g] Liquid mass in tank 3 at time t0
m40 = 9378.49308238605;                    % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 = 100;
F4_0 = 150;
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];           % [cm3/s] Disturbance variables at t0
% d_k = d0.*ones(2, length(t));
y0 = sensor_wo_noise(x0', At, rho);


% Modelling the disturbance as Brownian motion
Ns = length(d0); % Number of realizations
seed = 10;
[W,t,dW] = ScalarStdWienerProcess(tf,N,Ns,seed);
sigma = [2^2 0; 0 2^2];                             % Covariance for disturbances in F3 and F4
d_k = d0 + sigma*dW';

%% 

%linearization
% Steady State
xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);

[A, B, C, Cz, E, Gw] = system_matrices(xs, At, at, rho, gamma1, gamma2, g);
D = zeros(2,4);
%ZOH Discretization of Linear System
[Ad,Bd,Ed]=c2dzoh(A,B,E,dt);
sys = ss(Ad,Bd,C(1:2,:),D(:,1:2));

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;           % Covariance for process noise


% Augmentation after discritization
Ad_aug = [Ad, Ed; zeros(size(Ed, 2), size(Ad, 1)), eye(size(Ed, 2))];
Bd_aug = [Bd; zeros(size(Ed, 2), size(Bd, 2))];
Ed_aug = [Ed; zeros(size(Ed, 2), size(Ed, 2))];
C_aug = [C, zeros(size(C, 1), size(Ed, 2))];
Q_aug = [Q, zeros(size(Q, 1), size(Ed, 2)); zeros(size(Ed, 2), size(Q, 2)), eye(size(Ed, 2))];
Gw_aug =  [Gw, zeros(4,2); zeros(2,4) eye(2,2)];


% Design MPC
A = sys.A;
B = sys.B;
C = sys.C;

% Define MPC parameters
Qz = 10000 * eye(size(C, 1));  % Weight on output tracking
S = 1 * eye(size(B, 2));  % Weight on control effort

MPC_sys = UnconstrainedMPCDesign(A, B, C, Qz, S, Ph);

h_sp = [y0(1); y0(2)]; % [cm] Height set point
hdev_sp = h_sp - y0(1:2); % [cm] Height set point i deviation variable

% Set point at all times
% Rsp_dev = [hdev_sp(1)*ones(1, Ph); hdev_sp(2)*ones(1, Ph)];  % To referencer for Tank 1 og Tank 2
Rsp_dev = [hdev_sp(1) * ones(1,length(t)+Ph); 
       hdev_sp(2) * ones(1,length(t)+Ph)];


% Find tidspunkt for halvdelen af simulationen
stepchange1 = round(length(t)/4); 
stepchange2 = round(length(t)*2/4); 
% Opdater sætpunktet for tank 1 til 40 cm efter halvdelen af simulationen
Rsp_dev(1, stepchange1:end) = 60 - y0(1); 
Rsp_dev(2, stepchange2:end) = 70 - y0(2);

Rsp = [Rsp_dev(1,:) + y0(1); Rsp_dev(2,:) + y0(2)];

% d_k(1, stepchange1:end) = d_k(1, stepchange1:end) + 100; 
% d_k(2, stepchange2:end) = d_k(2, stepchange2:end) - 50;


x(:,1) = x0;
u(:,1) = u0;
for i =1:tf/dt
    % Beregn nuværende indeks og prediktionsindeks
    pred_idx = i:min(i+Ph-1, length(Rsp_dev)); % Sørg for ikke at gå ud over tidsvektoren
    % Skab prediktionshorisont for sætpunkt
    Rsp_dev_pred = Rsp_dev(:, pred_idx);   
    % Reshape til en 2Phx1 vektor
    ref_traj = reshape(Rsp_dev_pred, [], 1);

    xdev(:,i) = x(:,i) - x0;
    udev(:,i) = u(:,i) - u0;
    ddev(:,i) = d_k(:,i) - d0;

    % Kalman filter
    %--------------------------------------------------------------------------
    [x_hat, x_phat] = kalman_filter_aug_dynamic_pred(t(i), xdev(:,i), udev(:,i), ddev(:,i), At, rho, R, Q_aug, Ad_aug, Bd_aug, Ed_aug, Gw_aug, C_aug, Ph);
    %--------------------------------------------------------------------------
    x_mpc = [x_hat(1:4,1) x_phat(1:4,:)];

    for j = 1:Ph+1
        % Beregn den lineære del af omkostningsfunktionen
        g = MPC_sys.M_x0 * x_mpc(:,j) + MPC_sys.M_r * ref_traj;        
        % Løs QP-problemet for denne iteration
        u_current = qpsolver(MPC_sys.H, g, [], [], [], [], [], []);
        % Gem det første input fra løsningen
        u_mpc((j-1)*size(B, 2)+1:j*size(B, 2)) = u_current(1:size(B, 2));
    end
    % Konverter u_mpc til u_pred
    u_pred = reshape(u_mpc, 2, Ph+1) + u0;
 

    [T, X, D, U, x_discrete] = discrete_fourtankProcess_plus_noise(x(:,i), [t(i) t(i+1)], u(:,i), d_k(:,i), p, Q);
    x(:,i+1) = x_discrete(end,:)';


    % Apply first control action
    udev(:, i+1) = u_mpc(1:2); % First control input for the time step
    % udev(:, i+1) = u_mpc(end-1:end);
    u(:,i+1) = udev(:,i+1) + u0;
end

y = sensor_plus_noise(x',At,rho,R);



figure(1)
% Øvre subplot - Tankhøjde
for i = 1:2
    subplot(2,2,i)
    plot(t/60, y(i,:),'b', 'LineWidth', 2);
    hold on
    plot(t/60, Rsp(i,1:end-Ph),'k', 'LineWidth', 1);
    grid on;
    ylabel('height [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    ylim([40 90]);
    legend('Tank level', 'Set point', 'Location', 'best');
    title(['Tank ', num2str(i)], 'FontSize', 10);
end

% Nedre subplot - Flow og Setpunkt
for i = 1:2
    subplot(2,2,i+2)    
    % Primær y-akse for flow
    yyaxis left
    plot(t/60, u(i,:),'b', 'LineWidth', 2); 
    ylabel('[cm^3/s]', 'FontSize', 12);
    ylim([0 800]);    
    % Sekundær y-akse for sætpunkt
    yyaxis right
    plot(t/60, Rsp(i,1:end-Ph),'k', 'LineWidth', 1); 
    ylabel('Set point [cm]', 'FontSize', 12);
    ylim([40 90]);
    grid on;
    xlim([0 t(end)/60]);
    legend('Manipulated variable', 'Set point', 'Location', 'best');
    title(['F', num2str(i)], 'FontSize', 10);
end
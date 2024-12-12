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
tf= 60;                  % [s] End time
dt = 1;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
Ph = 40;                     % Prediction horizon

m10 = 17612.0123864868;                    % [g] Liquid mass in tank 1 at time t0
m20 = 29640.6694933624;                    % [g] Liquid mass in tank 2 at time t0
m30 = 4644.21948249842;                    % [g] Liquid mass in tank 3 at time t0
m40 = 9378.49308238599;                    % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 = 100;
F4_0 = 150;
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];           % [cm3/s] Disturbance variables at t0
d = d0.*ones(2, length(t));
u = u0.*ones(2, length(t));
[y0] = sensor_wo_noise(x0', At, rho);

%%  -------------------- 8.1 MPC function ------------------------------

%linearization
% Steady State
xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);
ys = sensor_wo_noise(xs,at,rho);
zs = sensor_wo_noise(xs,at,rho);

% %Stochastic Brownian
% sigma = [2^2 0; 0 2^2]; 
% [A, B, C, E, Gw] = linearized_models(xs, At, at, rho, gamma1, gamma2, g, ...
                                                                                % 'brownian', d0, sigma);

% --------------------------------------------------------------
% Linearization
% --------------------------------------------------------------
ys = mass_to_height(xs,At,rho);
hs = ys;
Tl = (At./at).*sqrt(2*hs/g);
sys.A=[-1/Tl(1) 0 1/Tl(3) 0;0 -1/Tl(2) 0 1/Tl(4);0 0 -1/Tl(3) 0;0 0 0 -1/Tl(4)];
sys.B=[rho*gamma1 0;0 rho*gamma2; 0 rho*(1-gamma2); rho*(1-gamma1) 0];
sys.C=diag(1./(rho*At));
Cz=sys.C(1:2,:);
E = [0 0; 0 0; rho 0; 0 rho];
Gw = eye(4);

A = sys.A;
B = sys.B;
C = sys.C;

%ZOH Discretization of Linear System
%Stochastic Brownian
[Ad,Bd,Ed]=c2dzoh(A,B,E,dt);
D = zeros(2,4);
% sys = ss(Ad,[Bd,Gd],C(1:2,:),D);
sys = ss(Ad,Bd,C(1:2,:),D(:,1:2));

% %Markov parameters
% [x11_brownian, x12_brownian, x21_brownian, x22_brownian] = MarkovPara(Ad,Bd,C,D,N);

%generate constants (this is moved inside "UnconstrainedMPCDesign")
% phi = generate_phi(A, C, Ph);
% Gamma = generate_Gamma(A, B, C, Ph);

% Define MPC parameters
Qz = 100 * eye(size(C, 1));  % Weight on output tracking
S = 0.1 * eye(size(B, 2)); % Weight on control effort
% Ph is Prediction horizon

% Design MPC
MPC_sys = UnconstrainedMPCDesign(A, B, E, C, Qz, S, Ph);


% Kalman filter parameters
R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*0.000004;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*0.0000004;           % Covariance for process noise


% ZOH Discretization of Linear System
[Ad, Bd, Gwd] = c2dzoh(A, B, Gw, dt);

% Augmentering efter diskretisering
Ad_aug = [Ad, Ed; zeros(size(Ed, 2), size(Ad, 1)), eye(size(Ed, 2))];
Bd_aug = [Bd; zeros(size(Ed, 2), size(Bd, 2))];
Ed_aug = [Ed; zeros(size(Ed, 2), size(Ed, 2))];
C_aug = [C, zeros(size(C, 1), size(Ed, 2))];
Q_aug = [Q, zeros(size(Q, 1), size(Ed, 2)); zeros(size(Ed, 2), size(Q, 2)), eye(size(Ed, 2))];
Gw_aug =  [Gwd, zeros(4,2); zeros(2,4) eye(2,2)];

% Generate noisy measurements
y_meas = C * x0 + 0.1 * randn(size(C, 1), tf / dt); % Noisy measurement example

% % Pack inputs
% sigma_R = R_hat; % Use the measurement noise covariance here
u0 = zeros(size(B, 2), 1); % Initial control input
Rd = 5;
% d_k = mvnrnd(zeros(tf/dt,size(C,1)),eye(2)*Rd)';
% v_k = mvnrnd(zeros(tf/dt,size(C,1)),eye(2)*sigma_R)';

% Cholesky decomposition of measurement noise covariance R
Lr = chol(R, 'lower');                  % Decompose R into lower triangular matrix
v_k = Lr * (randn(4,length(t)));           % Generate measurement noise with covariance R

R1 = [0 ; 0];
R2 = [10 ; 10];
R3 = [20 ;-10];

Rsp=zeros(2*(Ph+tf/dt),1);
Rsp(1:2*(tf/dt/3),1)=kron(ones(tf/dt/3,1),R1);
Rsp(2*(tf/dt/3)+1:2*(2*tf/dt/3),1)=kron(ones(tf/dt/3,1),R2);
Rsp(2*(2*tf/dt/3)+1:2*(Ph+tf/dt),1)=kron(ones(Ph+tf/dt/3,1),R3);

% Rsp = [20*ones(1,tf/dt);25*ones(1,tf/dt)];
inputs = {x0, u0, Rsp, d, v_k};

% Simulation
[y, u] = MPC_Sim_Unconstrained(sys, MPC_sys, Q_aug, Rsp, tf, dt, inputs, Ph,t,At,rho,Ad_aug, Bd_aug, Ed_aug, Gw_aug, C_aug, R);

figure(1)
for i = 1:4
    subplot(2,2,i)
    plot(t(1:60)/60, y(i,:),'b', 'LineWidth', 2); 
    hold off;
    grid on;
    ylabel('height [cm]', 'FontSize', 12);
    xlim([0 t(end)/60]);
    legend('Measured height', 'Dynamic Kalman filter', 'Static Kalman filter', 'Location', 'best');
    title(['Tank ', num2str(i)], 'FontSize', 10);
end
sgtitle('Not-Augmented Linear Kalman filter', 'FontSize', 14, 'FontWeight', 'bold');
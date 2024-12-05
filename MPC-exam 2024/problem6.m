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
tf= 400*60;                  % [s] End time
dt = 10;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
m10 = 17612.0123864868;                    % [g] Liquid mass in tank 1 at time t0
m20 = 29640.6694933624;                    % [g] Liquid mass in tank 2 at time t0
m30 = 4644.21948249842;                    % [g] Liquid mass in tank 3 at time t0
m40 = 9378.49308238599; 
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 =100;
F4_0 =150;
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];           % [cm3/s] Disturbance variables at t0
d = d0.*ones(2, length(t));
u = u0.*ones(2, length(t));


xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);    % Løser differentialligningssystemet
% (QuadrupleTankProcess) og solver hvad x er når hældningen er 0. Dvs. at
% den beregner hvad masserne er når der opnås steady state i tankene.

%%
% --------------------------------------------------------------
% Linearization
% --------------------------------------------------------------
ys = mass_to_height(xs,At,rho);
zs = output(xs',At,rho);
hs = ys;
T = (At./at).*sqrt(2*hs/g);
A=[-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B=[rho*gamma1 0;0 rho*gamma2; 0 rho*(1-gamma2); rho*(1-gamma1) 0];
C=diag(1./(rho*At));
Cz=C(1:2,:);


R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;     % Covariance for process noise


% Ns = length(d0); % Number of realizations
% seed = 100;
% [W,t,dW] = ScalarStdWienerProcess(tf,N,Ns,seed);
% 
% sigma = [2^2 0; 0 2^2];                             % Covariance for disturbances in F3 and F4
% d_brownian = d0 + sigma*dW';  

d(1,200:end) = d(1,200)*1.5;
d(2,800:end) = d(2,800)*0.75;

[T, X, D, U, x_sample] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y_sample] = sensor_plus_noise(x_sample, At, rho, R);


figure(1);
for i = 1:4
    subplot(2,2,i)
    plot(t/60, y_sample(i,:), 'LineWidth', 2);
    grid on;
    xlabel('t [min]', 'FontSize', 12);
    ylabel('h [cm]', 'FontSize', 12);
    xlim([0 T(end)/60]);
    title(['Liquid height in tank', num2str(i)], 'FontSize', 10);
end

[x_hat] = kalman_filter(t, x0, y_sample, At, rho, R, A, B, C)

% xk_k1 = x0; % Initial state
% Pk_k1 = 1000* eye(4); % Initial covariance
% x_hat = zeros(4, length(t)); % Pre-allocate estimated states
% y_hat = zeros(4, length(t)); % Pre-allocate estimated outputs
% x_pred_total =[];
% Lr = chol(R,'lower'); 
% for k = 1:length(t)
% 
%     % only filtering
%     yk_k1 = mass_to_height(xk_k1,At,rho);
%     ek = y_sample(:,k)-yk_k1;
%     Re_k=C*Pk_k1*C'+Lr;
%     K = Pk_k1*C'*Re_k^-1;
%     xk_k=xk_k1+K*ek;
%     Pk_k=Pk_k1-K*Re_k*K';
% 
%     x_hat(:,k) = xk_k1;
% 
%     xk_k1 =xk_k;
%     Pk_k1 = Pk_k;  
% end


figure(2)
for i = 1:4
    subplot(2,2,i)
    plot(t/60, x_sample(:,i),t/60, x_hat(i,:), 'LineWidth', 2);
    grid on;
    xlabel('t [min]', 'FontSize', 12);
    ylabel('m [g]', 'FontSize', 12);
    xlim([0 T(end)/60]);
    legend('Mass in tank', 'Kalman filter', 'Location', 'best');
    title(['Liquid height in tank', num2str(i)], 'FontSize', 10);
end
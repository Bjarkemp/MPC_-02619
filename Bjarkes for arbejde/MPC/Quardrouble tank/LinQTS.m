clc; clear all; close all;

%plan: løs systemet diskret og indsæt ss som gæt på lineariseringen

t0 = 0;
t_f = 1200; % 20 minutes in seconds
Ts = 4; % Sampling time [s]
t = t0:Ts:t_f; % Sampling instants [s]
num_steps = length(t);

% Initial parameters
m10 = 0.0; m20 = 0.0; m30 = 0.0; m40 = 0.0;
x0 = [m10; m20; m30; m40];
F1 = 250; % Initial F1
F2 = 325; % Constant F2

% Apply step increases at different times
u = [repmat(F1,1,num_steps); repmat(F2,1,num_steps)]; % Initial control input

% Define the time step at which the change happens (e.g., after 200 seconds)
step_time_index = find(t == 200); % Find the time index where the step happens

% Simulate 5%, 10%, and 25% step increases in F1
u(1, step_time_index:end) = F1 * 1.05;  % 5% increase in F1 from 200s onwards
u(1, step_time_index + 30:end) = F1 * 1.10;  % 10% increase after 250s
u(1, step_time_index + 60:end) = F1 * 1.25;  % 25% increase after 300s

umin = [0; 0];
umax = [400; 1000];

% Parameters
p = [1.2272; 1.2272; 1.2272; 1.2272; 380.1327; 380.1327; 380.1327; 380.1327; 981; 0.45; 0.40; 1];

% Pre-allocate for performance
X = zeros(num_steps, 4); % System states
T = zeros(num_steps, 1); % Time vector
U = zeros(num_steps, 2); % Control inputs (F1 and F2)

y = zeros(4, num_steps); % Measured states (with noise)
z = zeros(2, num_steps); % Simplified measured outputs (m1 and m2)

% Process Noise
Q = [20^2 0; 0 40^2]; % Process noise covariance matrix
Lq = chol(Q,'lower');
w = Lq * randn(2, num_steps); % Process noise for each time step

% Measurement Noise
R = eye(4); % Measurement noise covariance matrix
Lr = chol(R,'lower');
v = Lr * randn(4, num_steps); % Measurement noise for each time step

%Linearization

% Parameters
ap = p(1:4); % [cm2] Pipe cross sectional areas
At = p(5:8); % [cm2] Tank cross sectional areas
gam = p(10:11); % [-] Valve constants
g = p(9); %[cm/s2] The acceleration of gravity
rho = p(12); %[g/cm3] Density of water
p2 = [ap; At; gam; g; rho];
% Steady State
us = [F1; F2]; % [cm3/s] Flow rates
xs0 = [5000; 5000; 5000; 5000]; % [g] Initial guess on xs
xs = fsolve(@FourTankSystemWrap,xs0,[],us,p);
ys = FourTankSystemSensor(xs,p);
zs = FourTankSystemOutput(xs,p);


% Linearization

hs = ys;
T = (At./ap).*sqrt(2*hs/g);
A=[-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B=[rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
C=diag(1./(rho*At));
Cz=C(1:2,:);


% ZOH Discretization of Linear System

[Ad,Bd]=c2dzoh(A,B,Ts)

% Poles

p = eig(A)

%Transfer function
D = zeros(size(C, 1), size(B, 2));  % D is of size (number of outputs) x (number of inputs)
s = tf('s')
sys_ss = ss(A, B, C, D);  % Create state-space system
G = tf(sys_ss);  % Convert state-space to transfer function




% % Initial state
% x = x0;
% 
% for k = 1:num_steps-1
%     % Measurement with noise
%     y(:,k) = (x + v(:,k))/(p(12)*p(5)); % Measured levels with noise at time step k
%     z(:,k) = [x(1), x(2)]/(p(12)*p(5)); % Simplified measurement for first two tanks
% 
%     % Simulate process from t(k) to t(k+1)
%     [T_temp, X_temp] = ode15s(@(t,x) QuadrupleTankProcess(t, x, u(:,k) , p), [t(k) t(k+1)], x);
% 
%     % Update state
%     x = X_temp(end, :)'; % Final state after time step
%     T(k+1) = t(k+1); % Update time record
%     X(k+1, :) = x'; % Store state history   
% end
% 
% % Plot results
% figure;
% % Plot tank liquid mass (system response)
% subplot(2,1,1);
% plot(T, X, 'LineWidth', 1);
% xlabel('Time (s)');
% ylabel('Tank liquid mass (g)');
% legend('m1', 'm2', 'm3', 'm4');
% title('Quadruple Tank System Response (With Step Increases in F1)');
% xlim([0,1200]);
% ylim([0,16000]);
% 
% % Plot manipulated variables (control inputs)
% subplot(2,1,2);
% stairs(T, (u)', 'LineWidth', 1);
% xlabel('Time (s)');
% ylabel('Control Input (cm^3/s)');
% legend('F1', 'F2');
% title('Control Inputs over Time');
% xlim([0,1200]);
% ylim([0,500]);
% 
% % Plot results
% figure;
% % Plot tank liquid height (system response)
% subplot(2,2,1);
% plot(T, X(:,3)/(p(12)*p(7)), 'LineWidth', 1);
% xlabel('Time (s)');
% ylabel('Tank liquid height (cm)');
% legend('h3');
% title('Tank 3');
% xlim([0,1200]);
% ylim([0,50]);
% 
% % Plot tank liquid height (system response)
% subplot(2,2,2);
% plot(T, X(:,4)/(p(12)*p(8)), 'LineWidth', 1);
% xlabel('Time (s)');
% ylabel('Tank liquid height (cm)');
% legend('h4');
% title('Tank 4');
% xlim([0,1200]);
% ylim([0,50]);
% 
% % Plot tank liquid height (system response)
% subplot(2,2,3);
% plot(T, X(:,1)/(p(12)*p(5)), 'LineWidth', 1);
% xlabel('Time (s)');
% ylabel('Tank liquid height (cm)');
% legend('h1');
% title('Tank 1');
% xlim([0,1200]);
% ylim([0,50]);
% 
% % Plot tank liquid height (system response)
% subplot(2,2,4);
% plot(T, X(:,2)/(p(12)*p(6)), 'LineWidth', 1);
% xlabel('Time (s)');
% ylabel('Tank liquid height (cm)');
% legend('h2');
% title('Tank 2');
% xlim([0,1200]);
% ylim([0,50]);
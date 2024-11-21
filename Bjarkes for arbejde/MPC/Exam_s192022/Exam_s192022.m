clc; clear all; close all;

syms m1 m2 m3 m4 F1 F2 F3 F4 h1 h2 h3 h4 gamma1 gamma2 A1 A2 A3 A4 g rho a1 a2 a3 a4
%plan: løs systemet diskret og indsæt ss som gæt på lineariseringen

%% time
t0 = 0;
t_f = 1600; % 20 minutes in seconds
Ts = 4; % Sampling time [s]
t = t0:Ts:t_f; % Sampling instants [s]
num_steps = length(t);

%% Initial parameters
m10 = 0.0; m20 = 0.0; m30 = 0.0; m40 = 0.0; 
x0 = [m10; m20; m30; m40]; % tanks are empty when we start
F1 = 250; % Initial F1
F2 = 325; % initial F2
F3 = 200; % initial F3
F4 = 200; % initial F4

%% inputs
u = [repmat(F1,1,num_steps); repmat(F2,1,num_steps)]; % Initial control input
d_det = [repmat(F3,1,num_steps); repmat(F4,1,num_steps)]; % Disturbances

%% Generate Brownian Motion for Stochastic Disturbances using ScalarStdWienerProcess
T = t_f;     % Total time
N = 1000; % Number of time intervals
Ns = 2;      % Number of realizations (one for F3, one for F4)
seed = 123;  % Optional seed for reproducibility
sigma = 1000;

[W,Tw,dW] = ScalarStdWienerProcess(T, N, Ns, seed);

% Define the stochastic disturbances based on Wiener process
d_sto = [d_det(1, :) + sigma*W(1, num_steps); d_det(2, :) + sigma*W(2, num_steps)];  % Add initial disturbances to Wiener process
%% Bjarke du har rodet med det her, der er noget galt
%% conditions
umin = [0; 0];
umax = [400; 1000];

%% Parameters
p = [1.2272; 1.2272; 1.2272; 1.2272; 380.1327; 380.1327; 380.1327; 380.1327; 981; 0.45; 0.40; 1];

a1 = p(1); a2 = p(2); a3 = p(3); a4 = p(4);
    A1 = p(5); A2 = p(6); A3 = p(7); A4 = p(8);
    g = p(9); gamma1 = p(10); gamma2 = p(11); rho = p(12);
%% Linearization

% Parameters
ap = p(1:4); % [cm2] Pipe cross sectional areas
At = p(5:8); % [cm2] Tank cross sectional areas
gam = p(10:11); % [-] Valve constants
g = p(9); %[cm/s2] The acceleration of gravity
rho = p(12); %[g/cm3] Density of water
p2 = [ap; At; gam; g; rho];
% Steady State
us = [F1; F2];% [cm3/s] Flow rates
ds = [F3; F4];% [cm3/s] Flow rates
xs0 = [5000; 5000; 5000; 5000]; % [g] Initial guess on xs
xs = fsolve(@FourTankSystemWrap,xs0,[],us,ds,p);
ys = FourTankSystemSensor(xs,p);
zs = FourTankSystemOutput(xs,p);

hs = ys;
T = (At./ap).*sqrt(2*hs/g);
A=[-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B=[rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
C=diag(1./(rho*At));
Cz=C(1:2,:);

%% ZOH Discretization of Linear System

[Ad,Bd] = c2dzoh(A,B,Ts);

%% Transfer function
D = zeros(size(C, 1), size(B, 2));  % D is of size (number of outputs) x (number of inputs)
s = tf('s');
sys_ss = ss(A, B, C, D);  % Create state-space system
G = tf(sys_ss);  % Convert state-space to transfer function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Assignments %%%%%%%%%%%%%%%%%%%%%%%%%%%

%2.4

%% ODE Solver for Deterministic Model
% Pass deterministic disturbances (d_det)
[T_det, X_det] = ode15s(@(t,x) QuadrupleTankProcess(t,x,u(:,1),d_det(:,1),p), [t0 t_f], x0);

%% ODE Solver for Stochastic Model
% Simulate the system using stochastic disturbances (Brownian motion)
[T_sto, X_sto] = ode15s(@(t,x) QuadrupleTankProcess(t,x,u(:,1),d_sto(:,1),p), [t0 t_f], x0);

% Initialize Y_det and Y_sto to match the size of X_det and X_sto
Y_det = zeros(length(T_det), 4);
Y_sto = zeros(length(T_sto), 4);

% Convert mass (X) to height (Y) for the deterministic model
for k = 1:length(X_det)
    Y_det(k, :) = FourTankSystemSensor(X_det(k, :)', p)'; % Transpose to fill the row
end
% Convert mass (X) to height (Y) for the stochastic model
for k = 1:length(X_sto)
    Y_sto(k, :) = FourTankSystemSensor(X_sto(k, :)', p)'; % Transpose to fill the row
end

% Plotting Results

% Plot F1 and F2 in a loop
figure;
for i = 1:2
    subplot(2,1,i);
    plot(t, u(i,:), 'b');
    xlabel('Time (s)');
    ylabel('Flow rate (cm^3/s)');
    title(['Control Input F' num2str(i)]);
    xlim([0 t_f]);
end

% Plot deterministic and stochastic outputs for each tank in a loop
figure;
for i = 1:4
    subplot(2,2,i);
    plot(T_det, Y_det(:,i), 'b', T_sto, Y_sto(:,i), '--r');
    xlabel('Time (s)');
    ylabel('Height (cm)');
    legend('Deterministic', 'Stochastic');
    xlim([0 t_f]);
end


%% 4 - Step Responses
%% 4.1 - Step Responses for Deterministic Model with Specific Step Times
% Apply step changes to u using a loop

% Define the time intervals and corresponding step percentages
t_intervals = [1, 400/Ts, 800/Ts, 1200/Ts, num_steps]; % Time indices for each segment
percentages = [1, 1.10, 1.25, 1.50]; % Corresponding percentage increases

% Initialize arrays to store results
T_det_all = [];
X_det_all = [];
u_step = zeros(2, num_steps); % To track F1 and F2 inputs over time

% Initial conditions
x_current = x0;  % Start from the initial state
t_start = t0;

% Loop over each interval and apply the step changes
for i = 1:length(percentages)
    % Define the end time for the current interval
    t_end = t_intervals(i+1) * Ts;  % End time for this segment
    
    % Set the constant input for this interval
    u_current = [F1 * percentages(i); F2 * percentages(i)];
    
    % Solve ODE for the current segment with constant u_current
    [T_det, X_det] = ode15s(@(t,x) QuadrupleTankProcess(t, x, u_current, d_det(:,1), p), [t_start t_end], x_current);
    
    % Append results for time and state
    T_det_all = [T_det_all; T_det];
    X_det_all = [X_det_all; X_det];
    
    % Update u_step across the time interval for plotting inputs
    u_step(:, t_intervals(i):t_intervals(i+1)) = repmat(u_current, 1, t_intervals(i+1) - t_intervals(i) + 1);
    
    % Update initial conditions and start time for the next interval
    x_current = X_det(end, :);
    t_start = T_det(end);
end

% Convert mass (X) to height (Y) for the deterministic model
Y_det = zeros(length(T_det_all), 4);
for k = 1:length(X_det_all)
    Y_det(k, :) = FourTankSystemSensor(X_det_all(k, :)', p)'; % Transpose to fill the row
end

% Plotting Results

% Plot the control inputs F1 and F2 over time
figure;
for i = 1:2
    subplot(2,1,i);
    plot(t, u_step(i,:), 'b');
    xlabel('Time (s)');
    ylabel('Flow rate (cm^3/s)');
    legend(['F' num2str(i)],location='southeast');
    title(['Control Input F' num2str(i)]);
    xlim([0 t_f]);
end

% Plot deterministic outputs for each tank in a loop
figure;
for i = 1:4
    subplot(2,2,i);
    plot(T_det_all, Y_det(:,i), 'b');
    xlabel('Time (s)');
    ylabel('Height (cm)');
    title(['Height in Tank ' num2str(i)]);
    legend('Deterministic',location='southeast');
    xlim([0 t_f]);
end

%% 4.2 - Step Responses with noise
%make loppe with low, mid and high noise
% Process Noise
Q = [20^2 0 0 0; 0 20^2 0 0; 0 0 20^2 0; 0 0 0 20^2];
Lq = chol(Q,'lower'); 
w = Lq * randn(4, num_steps); % Process noise for each time step

% Measurement Noise
R = eye(4); 
Lr = chol(R,'lower');
v = Lr * randn(4, num_steps); % Measurement noise for each time step

% Allocate Arrays to Store Results
T_det_all = []; % Time
X_det_all = []; % States
u_step = zeros(2, num_steps); % Inputs for plotting

% Initialize State and Time for Loop
x_current = x0;
t_start = t0;

% Loop over each interval and apply the step changes
for i = 1:length(percentages)
    % Define the end time for the current interval
    t_end = t_intervals(i+1) * Ts;  % End time for this segment
    
    % Set the constant input for this interval
    u_current = [F1 * percentages(i); F2 * percentages(i)];
    
    % Solve ODE for the current segment with constant u_current
    [T_det, X_det] = ode15s(@(t,x) QuadrupleTankProcess(t, x, u_current, [0;0], p), [t_start t_end], x_current);
    
    % Add process noise to the states outside the ODE solver
    X_det = X_det + w(:,1:size(X_det,1))'; % Apply noise to each state
    
    % Append results for time and state
    T_det_all = [T_det_all; T_det];
    X_det_all = [X_det_all; X_det];
    
    % Update u_step across the time interval for plotting inputs
    u_step(:, t_intervals(i):t_intervals(i+1)) = repmat(u_current, 1, t_intervals(i+1) - t_intervals(i) + 1);

    % Update initial conditions and start time for the next interval
    x_current = X_det(end, :);
    t_start = T_det(end);
end

% Convert mass to height with measurement noise
Y_det = zeros(size(X_det_all,1), 4);
for k = 1:size(X_det_all,1)
    Y_det(k, :) = FourTankSystemSensor(X_det_all(k, :)', p)' + v(:,k)'; % Apply measurement noise to each output
end

% Plotting Results

% Plot the control inputs F1 and F2 over time
figure;
for i = 1:2
    subplot(2,1,i);
    plot(t, u_step(i,:), 'b');
    xlabel('Time (s)');
    ylabel('Flow rate (cm^3/s)');
    legend(['F' num2str(i)], location = 'southeast');
    title(['Control Input F' num2str(i)]);
    xlim([0 t_f]);
end

% Plot deterministic outputs for each tank in a loop
figure;
for i = 1:4
    subplot(2,2,i);
    plot(T_det_all, Y_det(:,i), 'b');
    xlabel('Time (s)');
    ylabel('Height (cm)');
    title(['Height in Tank ' num2str(i)]);
    legend('Deterministic with Noise', location = 'southeast');
    xlim([0 t_f]);
end


%% 4.3 Normalized step
% 
% % Storage for results
% step_u1 = {}; step_u2 = {};
% x_s = [5000; 5000; 5000; 5000]; % Assume some steady-state values
% 
% for j = 1:2 % Loop over both inputs F1 and F2
%     for i = 1:length(percentages)
%         % Apply step to F1 or F2 based on the outer loop `j`
%         if j == 1
%             u_current = [F1 * percentages(i); F2];
%         else
%             u_current = [F1; F2 * percentages(i)];
%         end
% 
%         % Solve the system with this step input
%         [T, X] = ode15s(@(t,x) QuadrupleTankProcess(t, x, u_current, [0;0], p), [t0:Ts:t_f], x0);
% 
%         % Compute normalized responses for h1 and h2
%         if j == 1 % Step in F1
%             Step_y1 = (X(:,1) - x_s(1)) / (u_current(1) - F1);
%             Step_y2 = (X(:,2) - x_s(2)) / (u_current(1) - F1);
%             step_u1{i} = {T, Step_y1, Step_y2};
%         else % Step in F2
%             Step_y1 = (X(:,1) - x_s(1)) / (u_current(2) - F2);
%             Step_y2 = (X(:,2) - x_s(2)) / (u_current(2) - F2);
%             step_u2{i} = {T, Step_y1, Step_y2};
%         end
%     end
% end
% 
% % Plot the normalized responses for F1 and F2 steps
% figure;
% % Plot normalized response for steps in F1
% subplot(2,2,1); hold on;
% for i = 1:length(percentages)
%     plot(step_u1{i}{1}/60, step_u1{i}{2}/(rho*A1), 'DisplayName', [num2str((percentages(i)-1)*100) '%']);
% end
% legend; hold off; grid on;
% xlabel('Time [min]'); ylabel('$h_1$ [Normalized]', 'Interpreter', 'latex');
% title('$u_1$ Step Response for $h_1$', 'Interpreter', 'latex');
% 
% subplot(2,2,3); hold on;
% for i = 1:length(percentages)
%     plot(step_u1{i}{1}/60, step_u1{i}{3}/(rho*A1), 'DisplayName', [num2str((percentages(i)-1)*100) '%']);
% end
% legend; hold off; grid on;
% xlabel('Time [min]'); ylabel('$h_2$ [Normalized]', 'Interpreter', 'latex');
% title('$u_1$ Step Response for $h_2$', 'Interpreter', 'latex');
% 
% % Plot normalized response for steps in F2
% subplot(2,2,2); hold on;
% for i = 1:length(percentages)
%     plot(step_u2{i}{1}/60, step_u2{i}{2}/(rho*A1), 'DisplayName', [num2str((percentages(i)-1)*100) '%']);
% end
% legend; hold off; grid on;
% xlabel('Time [min]'); ylabel('$h_1$ [Normalized]', 'Interpreter', 'latex');
% title('$u_2$ Step Response for $h_1$', 'Interpreter', 'latex');
% 
% subplot(2,2,4); hold on;
% for i = 1:length(percentages)
%     plot(step_u2{i}{1}/60, step_u2{i}{3}/(rho*A1), 'DisplayName', [num2str((percentages(i)-1)*100) '%']);
% end
% legend; hold off; grid on;
% xlabel('Time [min]'); ylabel('$h_2$ [Normalized]', 'Interpreter', 'latex');
% title('$u_2$ Step Response for $h_2$', 'Interpreter', 'latex');
% 
% %% 4.4 
% 
% % The steady-state gain for each transfer function can be calculated by 
% % observing the final value of the normalized response for each input-output pair. 
% 
% gain_u1y1 = step_u1{1}{2}(end) / (rho * A1);
% gain_u2y1 = step_u2{1}{2}(end) / (rho * A1);
% gain_u1y2 = step_u1{1}{3}(end) / (rho * A2);
% gain_u2y2 = step_u2{1}{3}(end) / (rho * A2);
% 
%% kalman filter

E = diag([0, 0, 2, 2]); % Each state has independent process noise; % Define the disturbance input matrix for process noise

% Define noise covariance matrices
Q = [20^2 0 0 0; 0 20^2 0 0; 0 0 20^2 0; 0 0 0 20^2]; % Process noise covariance
R = eye(4); % Measurement noise covariance, adjust if needed

% Initial State and Covariance
x_hat = [0; 0; 0; 0]; % Initial state estimate, change to see if filter is working
P = eye(4)*1000; % Initial error covariance estimate

% Store results for plotting
x_hat_all = zeros(4, num_steps); % Estimated states
y_all = zeros(4, num_steps); % Measured states
y_noisy_all = zeros(4, num_steps); % Measured outputs with noise

% Process and Measurement Noise
Lq = chol(Q, 'lower');
process_noise = Lq * randn(4, num_steps); % Process noise

Lr = chol(R, 'lower');
measurement_noise = Lr * randn(4, num_steps); % Measurement noise

% System Inputs
percentages = [1.1, 1.25, 1.5]; % Step percentages for F1 and F2
u_step = [repmat(250, 1, num_steps); repmat(325, 1, num_steps)]; % Adjust F1 and F2 as per steps

% Dynamic Kalman Filter Loop
for k = 1:num_steps-1
    % System Dynamics (Simulate actual system with process noise)
    u = u_step(:, k); % Inputs at time step k
    x_true = A * x_hat + B * u + E * process_noise(:, k); % True state update with process noise
    
    % Measurement Update (Simulated measurement with noise)
    y = C * x_true + measurement_noise(:, k);
    
    % Kalman Prediction Step
    x_hat_pred = A * x_hat + B * u; % Predicted state estimate
    P_pred = A * P * A' + Q; % Predicted error covariance
    
    % Kalman Gain Calculation
    K = P_pred * C' / (C * P_pred * C' + R); % Kalman gain
    
    % Kalman Update Step
    x_hat = x_hat_pred + K * (y - C * x_hat_pred); % Updated state estimate
    P = (eye(4) - K * C) * P_pred; % Updated error covariance
    
    % Store results
    x_hat_all(:, k) = x_hat; % Estimated state at time step k
    y_all(:, k) = C * x_hat; % Estimated output (without measurement noise)
    y_noisy_all(:, k) = y; % Noisy measurement
    
    % Update time step for next iteration
    x_hat = x_hat; % Update the state estimate for the next iteration
end

% Plotting Results
% Compare Kalman Filter estimates with true and noisy measurements

% Plot Estimated vs Noisy Measurements for each Tank
figure;
for i = 1:4
    subplot(2,2,i);
    plot(t, y_noisy_all(i,:), 'r--', 'DisplayName', 'Noisy Measurement');
    hold on;
    plot(t, y_all(i,:), 'b-', 'DisplayName', 'Kalman Estimate');
    xlabel('Time (s)');
    ylabel(['Height in Tank ' num2str(i) ' (cm)']);
    title(['Kalman Filter Estimate vs Measurement - Tank ' num2str(i)]);
    legend('Location', 'best');
    hold off;
    xlim([0 t_f]);
end

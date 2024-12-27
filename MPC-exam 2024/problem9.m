clc, clear, close all
addpath("Functions"); % Add the "Functions" folder to the MATLAB path

% Problem 8
% ----------------------------------------------------------
% Parameters
% ----------------------------------------------------------
a1 = 1.2272;                % [cm^2] Outlet pipe area for Tank 1
a2 = 1.2272;                % [cm^2] Outlet pipe area for Tank 2
a3 = 1.2272;                % [cm^2] Outlet pipe area for Tank 3
a4 = 1.2272;                % [cm^2] Outlet pipe area for Tank 4
A1 = 380.1327;              % [cm^2] Cross-sectional area for Tank 1
A2 = 380.1327;              % [cm^2] Cross-sectional area for Tank 2
A3 = 380.1327;              % [cm^2] Cross-sectional area for Tank 3
A4 = 380.1327;              % [cm^2] Cross-sectional area for Tank 4
g = 981;                    % [cm/s^2] Gravitational acceleration
gamma1 = 0.6;               % Flow distribution coefficient for valve 1
gamma2 = 0.7;               % Flow distribution coefficient for valve 2
rho = 1.0;                  % [g/cm^3] Density of water

% Parameter vector containing all constants
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho]; 
at = p(1:4);                % Outlet pipe areas [cm^2]
At = p(5:8);                % Cross-sectional areas [cm^2]

% -----------------------------------------------------------
% Simulation scenario
% -----------------------------------------------------------
t0 = 0.0;                   % [s] Initial time
tf = 60*20;                 % [s] End time (20 minutes)
dt = 1;                    % [s] Time step size
N = tf/dt;                  % Total number of steps
t = t0:dt:tf;               % [s] Time vector
Ph = 40;                    % Prediction horizon for MPC  

% Initial flow rates [cm^3/s]
F1_0 = 300;                 % Pump flow rate 1
F2_0 = 300;                 % Pump flow rate 2
F3_0 = 100;                 % External disturbance flow 1
F4_0 = 150;                 % External disturbance flow 2

u0 = [F1_0; F2_0];          % Initial manipulated variables (pumps)
d0 = [F3_0; F4_0];          % Initial disturbances
% d_k = d0.*ones(2,length(t));
% Initial tank masses [g]
xs0 = [5000; 5000; 5000; 5000]; % [g] Initial guess on xs
x0 = fsolve(@FourTankSystemWrap,xs0,[],u0,d0,p);

% Simulated sensor measurements without noise
y0 = sensor_wo_noise(x0', At, rho);

% -----------------------------------------------------------
% Modeling disturbance as Brownian motion
% -----------------------------------------------------------
Ns = length(d0);            % Number of realizations for disturbances
seed = 10;                  % Random seed for repeatability
[W, t, dW] = ScalarStdWienerProcess(tf, N, Ns, seed); % Generate Wiener process
sigma = [2^2 0; 0 2^2];     % Covariance matrix for disturbances
d_k = d0 + sigma * dW';     % Final disturbance vector

%% -----------------------------------------------------------
% Linearization at steady state
% -----------------------------------------------------------
xs = fsolve(@FourTankSystemWrap, x0, [], u0, d0, p); % Solve for steady state

% System matrices (linearized system)
[A, B, C, Cz, E, Gw] = system_matrices(xs, At, at, rho, gamma1, gamma2, g);
D = zeros(2,4);             % Direct feedthrough matrix (zeros)

% Discretization using Zero-Order Hold (ZOH)
[Ad, Bd, Ed] = c2dzoh(A, B, E, dt);

% Covariance matrices for process and measurement noise
R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4; % Measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;       % Process noise

% -----------------------------------------------------------
% System augmentation for Kalman filter
% -----------------------------------------------------------
Ad_aug = [Ad, Ed; zeros(size(Ed, 2), size(Ad, 1)), eye(size(Ed, 2))];
Bd_aug = [Bd; zeros(size(Ed, 2), size(Bd, 2))];
Ed_aug = [Ed; zeros(size(Ed, 2), size(Ed, 2))];
C_aug = [C, zeros(size(C, 1), size(Ed, 2))];
Q_aug = [Q, zeros(size(Q, 1), size(Ed, 2)); zeros(size(Ed, 2), size(Q, 2)), eye(size(Ed, 2))];
Gw_aug = [Gw, zeros(4,2); zeros(2,4), eye(2)];

% -----------------------------------------------------------
% Design MPC
% -----------------------------------------------------------
sys = ss(Ad, Bd, C(1:2,:), D(:,1:2)); % Discrete-time state-space system
Qz = 300 * eye(size(sys.C, 1));% Weight for output tracking
S = 1 * eye(size(sys.B, 2));    % Weight for control effort

%Constraints umin =< u =< umax !IN DEVIATION VARIABLES!
umin = 0; 
umax = 500;
Umin = repmat(umin - u0,2*Ph,1);
Umax = repmat(umax - u0,2*Ph,1);
lb_rate = repmat(-20,2*Ph,1);  % minimum rate of change
ub_rate = repmat(20,2*Ph,1);   % maximum rate  of change


%setting up system
MPC_sys = UnconstrainedMPCDesign(sys.A, sys.B, sys.C, Qz, S, Ph);

% -----------------------------------------------------------
% Setpoint trajectories with step changes
% -----------------------------------------------------------
h_sp = [y0(1); y0(2)];                % Initial setpoints for tank heights
hdev_sp = h_sp - y0(1:2);             % Setpoints in deviation form

Rsp_dev = [hdev_sp(1) * ones(1, length(t)+Ph); 
           hdev_sp(2) * ones(1, length(t)+Ph)];

% Step changes at specific time steps
stepchange1 = round(length(t)/4);     % First setpoint step change
stepchange2 = round(length(t)*2/4);   % Second setpoint step change
stepchange3 = round(length(t)*3/4);   % Third setpoint step change
Rsp_dev(1, stepchange1:end) = 60 - y0(1); % Updated setpoint for Tank 1
% Rsp_dev(2, stepchange2:end) = 70 - y0(2); % Updated setpoint for Tank 2
Rsp_dev(:, stepchange3:end) = 0; % Updated setpoint for Tank 2

Rsp = [Rsp_dev(1,:) + y0(1); Rsp_dev(2,:) + y0(2)]; % Final setpoint trajectories

% d_k(1, stepchange1:end) = d_k(1, stepchange1:end) + 100; 
d_k(2, stepchange2:end) = d_k(2, stepchange2:end) - 100;



% -----------------------------------------------------------
% Simulation Loop
% -----------------------------------------------------------
x(:,1) = x0; % Initialize system states
u(:,1) = u0; % Initialize control inputs;

for i = 1:tf/dt
   
    pred_idx = i:min(i+Ph-1, length(Rsp_dev)); % Prediction indices
    Rsp_dev_pred = Rsp_dev(:, pred_idx);       % Predicted setpoints
    ref_traj = reshape(Rsp_dev_pred, [], 1);   % Reshape into a single vector

    xdev(:,i) = x(:,i) - x0;   % State deviation
    udev(:,i) = u(:,i) - u0;   % Input deviation
    ddev(:,i) = d_k(:,i) - d0; % Disturbance deviation

    % Kalman filter for state estimation (Det kan godt være at vi ikke
    % behøver et filter der predicter, siden jeg tror der også sker noget
    % prediction i g
    [x_hat, x_phat] = kalman_filter_aug_dynamic_pred(t(i), xdev(:,i), udev(:,i), At, rho, R, Q_aug, Ad_aug, Bd_aug, Gw_aug, C_aug, Ph);
    x_mpc = [x_hat(1:4,1) x_phat(1:4,:)]; % Combine estimated states
    if i == 1 % can have a zero'th iteration on MV's
        dub = ub_rate + MPC_sys.I0 * udev(:,i);
        dlb = lb_rate + MPC_sys.I0 * udev(:,i);
        g = MPC_sys.M_x0 * x_hat(1:4,1) + MPC_sys.M_r * ref_traj + MPC_sys.M_u1 * udev(:,i); % Linear cost function
        u_current = qpsolver(MPC_sys.H, g, Umin, Umax, MPC_sys.Lambda, dlb, dub, x_hat(1:4,1)); % Solve QP problem
        u_pred = reshape(u_current', 2, Ph) + u0; % Predicted control inputs
        udev(:, i+1) = u_pred(:,1) - u0; % First control input for the time step
    else
        dub = ub_rate + MPC_sys.I0 * udev(:,i-1);
        dlb = lb_rate + MPC_sys.I0 * udev(:,i-1);
        g = MPC_sys.M_x0 * x_hat(1:4,1) + MPC_sys.M_r * ref_traj + MPC_sys.M_u1 * udev(:,i-1); % Linear cost function
        u_current = qpsolver(MPC_sys.H, g, Umin, Umax, MPC_sys.Lambda, dlb, dub, x_hat(1:4,1)); % Solve QP problem
        u_pred = reshape(u_current', 2, Ph) + u0; % Predicted control inputs
        udev(:, i+1) = u_pred(:,1) - u0; % First control input for the time step        
    end % if
    % u_pred = reshape(u_mpc, 2, Ph+1) + u0; % Predicted control inputs
    [~, ~, ~, ~, x_discrete] = discrete_fourtankProcess_plus_noise(x(:,i), [t(i) t(i+1)], u(:,i), d_k(:,i), p, Q);
    x(:,i+1) = x_discrete(end,:)'; % Update system state

    % Apply first control action
    % udev(:, i+1) = u_mpc(1:2); % First control input for the time step
    u(:,i+1) = udev(:,i+1) + u0;
end

% Simulate sensor measurements with added noise
y = sensor_plus_noise(x', At, rho, R);

% -----------------------------------------------------------
% Plot Results
% -----------------------------------------------------------

figure(1)

% Upper subplot - Tank levels
for i = 1:2
    subplot(2,2,i) % Create subplot for Tank i
    plot(t/60, y(i,:),'b', 'LineWidth', 2); % Plot measured tank levels
    hold on
    plot(t/60, Rsp(i,1:end-Ph),'k', 'LineWidth', 1); % Plot setpoint trajectory
    grid on; % Add grid to the plot
    ylabel('Height [cm]', 'FontSize', 12); % Y-axis label
    xlim([0 t(end)/60]); % Set x-axis limits (time in minutes)
    ylim([40 90]); % Set y-axis limits
    legend('Tank level', 'Set point', 'Location', 'best'); % Add legend
    title(['Tank ', num2str(i)], 'FontSize', 10); % Title for the subplot
end

% Lower subplot - Flow rate and setpoint
for i = 1:2
    subplot(2,2,i+2) % Create subplot for manipulated variables and setpoints

    % Primary Y-axis for flow rate
    yyaxis left
    plot(t/60, u(i,:),'b', 'LineWidth', 2); % Plot control inputs (flow rates)
    hold on;
    yline(umin, '--r', 'LineWidth', 1.5); % Plot umin as a dashed red line
    yline(umax, '--r', 'LineWidth', 1.5); % Plot umax as a dashed green line
    ylabel('[cm^3/s]', 'FontSize', 12); % Label for primary Y-axis (flow rate)
    ylim([0 800]); % Set limits for flow rate

    % Secondary Y-axis for setpoint
    yyaxis right
    plot(t/60, Rsp(i,1:end-Ph),'k', 'LineWidth', 1); % Plot setpoint trajectory
    ylabel('Set point [cm]', 'FontSize', 12); % Label for secondary Y-axis (setpoint)
    ylim([40 90]); % Set limits for setpoint

    grid on; % Add grid to the plot
    xlim([0 t(end)/60]); % Set x-axis limits (time in minutes)
    legend('Manipulated variable', 'umin', 'umax', 'Set point', 'Location', 'best'); % Add legend
    title(['F', num2str(i)], 'FontSize', 10); % Title for the subplot
end

%% save figure
saveas(gcf,fullfile('C:\Users\bjark\OneDrive\Skrivebord\MPC_-02619\MPC-exam 2024\Plots','problem9.tiff'),'tiff')
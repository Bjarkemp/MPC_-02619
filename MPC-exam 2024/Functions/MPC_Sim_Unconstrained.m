function [y, u, x_hat] = MPC_Sim_Unconstrained(A, B, C, MPC_sys, Q_hat, R_hat, tf, dt, inputs, pH)
    % MPC_SIM_UNCONSTRAINED simulates an MPC with a Kalman filter for state estimation.
    %
    % Inputs:
    %   A, B, C - Discrete-time state-space matrices
    %   MPC_sys - Struct containing MPC matrices (H, M_x0, M_r, etc.)
    %   Q_hat   - Process noise covariance matrix
    %   R_hat   - Measurement noise covariance matrix
    %   tf      - Final simulation time
    %   dt      - Sampling time
    %   inputs  - Cell array containing:
    %             {x0, u_vec, R, y_meas}
    %   pH      - Prediction horizon
    %
    % Outputs:
    %   y       - Output trajectory
    %   u       - Control input trajectory
    %   x_hat   - State estimate trajectory

    % Unpack MPC matrices
    H = MPC_sys.H;
    M_x0 = MPC_sys.M_x0;
    M_r = MPC_sys.M_r;

    % Unpack inputs
    x0 = inputs{1};    % Initial state
    u_vec = inputs{2}; % Initial control inputs (zeros)
    R = inputs{3};     % Reference trajectory
    y_meas = inputs{4}; % Measured outputs with noise

    % Simulation time steps
    num_steps = round(tf / dt);

    % Initialize variables
    x = zeros(size(A, 1), num_steps + 1);  % True states
    x_hat = zeros(size(A, 1), num_steps + 1); % State estimates
    u = zeros(size(B, 2), num_steps);      % Control inputs
    y = zeros(size(C, 1), num_steps);      % Outputs
    x(:, 1) = x0;                          % Set initial state
    x_hat(:, 1) = x0;                      % Initialize state estimate

    % Kalman filter initialization
    P = eye(size(A)); % Initial estimation error covariance

    % Ensure R is properly extended to cover the prediction horizon
    if size(R, 1) < pH * size(C, 1)
        R = repmat(R, ceil(pH * size(C, 1) / size(R, 1)), 1); % Extend if too short
    end

    % Loop over time steps
    for i = 1:num_steps
        % Compute reference for the current prediction horizon
        rows_remaining = min(size(R, 1) - (i-1)*size(C,1), pH*size(C,1));
        R_current = R((i-1)*size(C,1)+1:(i-1)*size(C,1)+rows_remaining, 1);

        % Solve unconstrained MPC problem
        g = M_x0 * x_hat(:, i) + M_r(1:size(H, 1), 1:rows_remaining) * R_current; % Linear term in the cost function
        u_mpc = -H \ g; % Solve quadratic program for unconstrained MPC (direct solve)

        % Apply first control input
        u(:, i) = u_mpc(1:size(B, 2));

        % Simulate system for the next state (true system dynamics)
        x(:, i+1) = A * x(:, i) + B * u(:, i);

        % Simulate true output
        y(:, i) = C * x(:, i);

        % Kalman filter: Measurement update
        K = P * C' / (C * P * C' + R_hat); % Kalman gain
        x_hat(:, i) = x_hat(:, i) + K * (y_meas(:, i) - C * x_hat(:, i)); % Update state estimate
        P = (eye(size(A)) - K * C) * P; % Update error covariance

        % Kalman filter: Time update
        x_hat(:, i+1) = A * x_hat(:, i) + B * u(:, i); % Predict next state
        P = A * P * A' + Q_hat; % Predict next error covariance
    end
end
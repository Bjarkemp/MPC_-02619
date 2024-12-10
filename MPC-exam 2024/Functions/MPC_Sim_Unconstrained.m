function [y, u] = MPC_Sim_Unconstrained(sys, MPC_sys, Q_hat, R_hat, t_f, Ts, inputs, Ph)
    % Inputs unpacking
    x = inputs{1};  % Initial state
    u = inputs{2};  % Control input (zeros initially)
    R = inputs{3};  % Reference trajectory
    d_k = inputs{4}; % Disturbances
    v_k = inputs{5}; % Measurement noise

    % System matrices
    A = sys.A; 
    B = sys.B; 
    E = sys.B(:, end-size(d_k, 1)+1:end); % Disturbance mapping
    C = sys.C;

    % Kalman filter initialization
    P = eye(size(A)); % Initial error covariance
    L = dlqe(A, eye(size(A)), C, Q_hat, R_hat); % Discrete-time Kalman gainR_hat

    % Initialize variables
    num_steps = t_f / Ts; % Total simulation steps
    x_hat = zeros(size(A, 1), num_steps + 1); % State estimate
    y = zeros(size(C, 1), num_steps);         % Output trajectory
    u = zeros(size(B, 2), num_steps);         % Control inputs
    

    % Loop over simulation
    for i = 1:num_steps
        % Extract reference trajectory for the prediction horizon
        rows_remaining = min(size(R, 1) - (i-1)*size(C,1), Ph*size(C,1));
        R_current = R((i-1)*size(C,1)+1 : (i-1)*size(C,1)+rows_remaining, 1);

        % If fewer rows than expected, pad with zeros
        if size(R_current, 1) < Ph * size(C, 1)
             R_current = [R_current; zeros(Ph * size(C, 1) - size(R_current, 1), 1)];
        end

        % Solve unconstrained MPC problem
        g = MPC_sys.M_x0 * x_hat(:, i) + MPC_sys.M_r * R_current;
        u_mpc = qpsolver(MPC_sys.H, g, [], [], [], [], [], []);

        % Apply first control action
        u(:, i) = u_mpc(1:size(B, 2)); % First control input for the time step

        % State and output updates
        x(:, i+1) = A * x(:, i) + B * u(:, i) + E * d_k(:, i); % State update
        y(:, i) = C * x(:, i) + v_k(:, i);                    % Output update

        % Kalman filter: Measurement update
        x_hat(:, i) = x_hat(:, i) + L * (y(:, i) - C * x_hat(:, i));

        % Kalman filter: Time update
        x_hat(:, i+1) = A * x_hat(:, i) + B * u(:, i);
    end
end

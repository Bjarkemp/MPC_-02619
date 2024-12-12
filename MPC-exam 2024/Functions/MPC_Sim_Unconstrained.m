function [y, u] = MPC_Sim_Unconstrained(sys, MPC_sys, Q_aug, Rsp, t_f, Ts, inputs, Ph,t,At,rho,Ad_aug, Bd_aug, Gd_aug, Gw_aug, C_aug,R)
    % Inputs unpacking
    x0 = inputs{1};  % Initial state
    x = x0;
    u0 = inputs{2};  % Control input (zeros initially)
    R = inputs{3};  % Reference trajectory
    d_k = inputs{4}; % Disturbances
    v_k = inputs{5}; % Measurement noise

    % System matrices
    A = sys.A; 
    B = sys.B; 
    E = sys.B(:, end-size(d_k, 1)+1:end); % Disturbance mapping
    C = sys.C;

    % % Kalman filter initialization
    % P = eye(size(A)); % Initial error covariance
    % L = dlqe(A, eye(size(A)), C, Q_hat, R_hat); % Discrete-time Kalman gainR_hat

    % Initialize variables
    num_steps = t_f / Ts; % Total simulation steps
    x_hat = zeros(size(A, 1), num_steps + 1); % State estimate
    y = zeros(size(C, 1), num_steps);         % Output trajectory
    u = zeros(size(B, 2), num_steps);         % Control inputs

    % Loop over simulation
    for i = 1:60%num_steps
        % Extract reference trajectory for the prediction horizon
        rows_remaining = min(size(Rsp, 1) - (i-1)*size(C,1), Ph*size(C,1));
        R_current = Rsp((i-1)*size(C,1)+1 : (i-1)*size(C,1)+rows_remaining, 1);

        % If fewer rows than expected, pad with zeros
        if size(R_current, 1) < Ph * size(C, 1)
             R_current = [R_current; zeros(Ph * size(C, 1) - size(R_current, 1), 1)];
        end
%--------------------------------------------------------------------------
        [x_hat2_dyn_pre, x_phat2_dyn_pre] = kalman_filter_aug_dynamic_pred(t(i), x0, u(:,i), d_k(:,i), At, rho, R, Q_aug, Ad_aug, Bd_aug, Gd_aug, Gw_aug, C_aug,Ph);
        x_hat = x_phat2_dyn_pre(1:4,:);
%--------------------------------------------------------------------------
        % Solve unconstrained MPC problem
        g = MPC_sys.M_x0 * x_hat(:, i) + MPC_sys.M_r * R_current;
        u_mpc = qpsolver(MPC_sys.H, g, [], [], [], [], [], []);

        % Apply first control action
        u(:, i) = u_mpc(1:size(B, 2)); % First control input for the time step

        % % State and output updates
        % x(:, i+1) = A * x(:, i) + B * u(:, i) + E * d_k(:, i); % State update
        y(:, i) = C * x_hat(:, i) + v_k(:, i);                     % Output update
        % 
        % % Kalman filter: Measurement update
        % x_hat(:, i) = x_hat(:, i) + L * (y(:, i) - C * x_hat(:, i));
        % 
        % % Kalman filter: Time update
        % x_hat(:, i+1) = A * x_hat(:, i) + B * u(:, i);



    end
end

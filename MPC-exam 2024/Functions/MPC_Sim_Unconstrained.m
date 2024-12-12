function [y, u] = MPC_Sim_Unconstrained(sys, MPC_sys, Q_aug, Rsp, tf, dt, inputs, Ph,t,At,rho,Ad_aug, Bd_aug, Gd_aug, Gw_aug, C_aug, R)
    % Inputs unpacking
    x0 = inputs{1};  % Initial state
    x = x0;
    u0 = inputs{2};  % Control input (zeros initially)
    % R = inputs{3};  % Reference trajectory
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
    num_steps = tf / dt; % Total simulation steps
    x_hat = zeros(size(A, 1), num_steps + 1); % State estimate
    y = zeros(size(C, 1), num_steps);         % Output trajectory
    u = zeros(size(B, 2), num_steps);         % Control inputs
    

    % Loop over simulation
    for i = 1:tf/dt
%--------------------------------------------------------------------------
        [x_hat2_dyn_pre, x_phat2_dyn_pre] = kalman_filter_aug_dynamic_pred(t(i), x0, u(:,i), d_k(:,i), At, rho, R, Q_aug, Ad_aug, Bd_aug, Gd_aug, Gw_aug, C_aug,Ph);
        x_hat = x_phat2_dyn_pre(1:4,:)
%--------------------------------------------------------------------------
        if i == 1
            u_mpc = Uncon_MPC(MPC_sys,x_hat(:,i),Rsp(1:2*Ph,1),u(:,i))
        else
            u_mpc = Uncon_MPC(MPC_sys,x_hat(:,i), Rsp(2*i+1:2*N+2*i,1), u(:,i-1))
        end
        u(:,i) = u_mpc(1:2,1); % Can only use from 1 to number of inputs

        % Apply first control action
        u(:, i) = u_mpc(1:size(B, 2)); % First control input for the time step
        y(:, i) = C * x_hat(:, i) + v_k(:, i);                     % Output update
    end
    u = u(:,1:end-1);
end

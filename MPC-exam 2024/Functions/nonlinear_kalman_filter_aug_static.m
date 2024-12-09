function [x_hat, x_phat] = nonlinear_kalman_filter_aug_static(t, xdev, udev, ddev, x0, u, d, At, rho, R, Q_aug, Ad_aug, Gw_aug, C_aug, p)
% nonlinear_kalman_filter_aug_static: Implements a nonlinear augmented Kalman filter
% with a static gain for a four-tank process system.

% Initialize variables
x_hat = [];
Lr = chol(R,'lower');                % Cholesky decomposition of measurement noise covariance
v = Lr * randn(size(xdev));          % Generate measurement noise
xhat_k_k1 = [xdev(:,1); ddev(:,1)];  % Initial augmented state estimate

% Solve the discrete-time algebraic Riccati equation for steady-state error covariance
P = idare(Ad_aug', C_aug', Gw_aug * Q_aug * Gw_aug', R); 

% Kalman filter loop
for k = 1:length(t)-1
    % ** Filtering Step **
    % Compute the innovation covariance
    Re_k = C_aug * P * C_aug' + R;

    % Compute Kalman gain
    K = P * C_aug' / Re_k;

    % Current measurement with noise
    yk = mass_to_height(xdev(:,k), At, rho) + v(:,k); 

    % Predicted measurement with noise
    yhat_k_k1 = mass_to_height(xhat_k_k1(1:4), At, rho) + v(:,k); 

    % Compute the innovation (measurement residual)
    ek = yk - yhat_k_k1; 

    % Update step: Correct state estimate using Kalman gain
    xhat_k_k = xhat_k_k1 + K * ek; 

    % ** Prediction Step **
    % Simulate the nonlinear system dynamics
    x_k_k = xhat_k_k(1:4) + x0;  % Corrected state in the original domain
    [~, Xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], x_k_k, [], u(:,k), d(:,k), p);

    % Compute deviation for next step
    x_phat(:,k) = [Xk(end,:)' - x0; d(:,k) - d(:,1)];

    % Store the corrected state
    x_hat = [x_hat xhat_k_k];  

    % Prepare for the next iteration
    xhat_k_k1 = x_phat(:,k);  
end

% Final update for the last time step
k = length(t);
Re_k = C_aug * P * C_aug' + R;
K = P * C_aug' / Re_k;
yk = mass_to_height(xdev(:,k), At, rho) + v(:,k); 
yhat_k_k1 = mass_to_height(xhat_k_k1(1:4), At, rho) + v(:,k); 
ek = yk - yhat_k_k1;  
xhat_k_k = xhat_k_k1 + K * ek; 

% Store the final corrected state
x_hat = [x_hat xhat_k_k];
end

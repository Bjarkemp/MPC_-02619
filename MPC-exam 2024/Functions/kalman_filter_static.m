function [x_hat, x_phat] = kalman_filter_static(t, xdev, udev, ddev, At, rho, R, Q, Ad, Bd, Gd, C)

% Compute Cholesky decomposition of the measurement noise covariance R
Lr = chol(R, 'lower');              % Decomposes R to lower triangular matrix
v = Lr * (randn(size(xdev)));       % Generate measurement noise with standard deviation given by R

% Initialize the state estimate with the first deviation value
xhat_k_k1 = xdev(:, 1);             % Initial state estimate

% Solve the algebraic Riccati equation to calculate steady-state covariance
P = idare(Ad', C', Q, R);           % Solve the Discrete Algebraic Riccati Equation for steady-state P

% Preallocate space for estimated states and predictions
x_hat = zeros(size(xdev));          % State estimates
x_phat = zeros(size(xdev));         % One-step-ahead state predictions

for k = 1:length(t)
    % Measurement prediction
    yhat_k_k1 = C * xhat_k_k1 ;   % Predicted measurement with noise
    ydev(:, k) = C * xdev(:, k) + v(:, k); % Actual measurement with noise

    % Innovation sequence (measurement residual)
    ek = ydev(:, k) - yhat_k_k1;         % Innovation between actual and predicted measurement

    % Innovation covariance
    Re_k = C * P * C' + R;               % Combine predicted error covariance and measurement noise

    % Kalman gain calculation
    K = P * C' / Re_k;                   % Optimal Kalman gain

    % Update step: update state estimate using the Kalman gain
    xhat_k_k = xhat_k_k1 + K * ek;       % Correct the prior estimate
    x_hat(:, k) = xhat_k_k;              % Store the corrected state estimate

    % One-step prediction: predict next state based on system dynamics
    x_phat(:, k) = Ad * xhat_k_k + Bd * udev(:, k) + Gd * ddev(:, k); % Predicted next state
    xhat_k_k1 = x_phat(:, k);           % Prepare the next iteration with predicted state
end

end

function [x_hat, x_phat] = kalman_filter_dynamic(t, xdev, udev, ddev, At, rho, R, Q, Ad, Bd, Gd, C)

% Compute Cholesky decomposition of measurement noise covariance R
Lr = chol(R, 'lower');              % Decompose R into a lower triangular matrix
v = Lr * (randn(size(xdev)));       % Generate measurement noise with covariance R

% Initialize the state estimate and error covariance
xhat_k_k1 = xdev(:, 1);             % Initial state estimate
P_k_k1 = 1000 * eye(size(Ad, 1));   % Initial error covariance, large to reflect high uncertainty

% Preallocate space for estimated states and predictions
x_hat = zeros(size(xdev));          % State estimates
x_phat = zeros(size(xdev));         % One-step-ahead state predictions

for k = 1:length(t)
    % Filtering step: predict and correct
    
    % Predict the measurement
    yhat_k_k1 = C * xhat_k_k1 + v(:, k);  % Predicted measurement
    ydev(:, k) = C * xdev(:, k) + v(:, k); % Actual measurement with noise

    % Compute the innovation sequence (measurement residual)
    ek = ydev(:, k) - yhat_k_k1;         % Difference between actual and predicted measurement

    % Compute the innovation covariance
    Re_k = C * P_k_k1 * C' + R;          % Combine predicted error covariance and measurement noise

    % Calculate the Kalman gain
    K = P_k_k1 * C' / Re_k;              % Optimal Kalman gain

    % Update step: correct the state estimate using the Kalman gain
    xhat_k_k = xhat_k_k1 + K * ek;       % Correct the prior estimate
    P_k_k = P_k_k1 - K * Re_k * K';      % Update the error covariance
    x_hat(:, k) = xhat_k_k;              % Store the corrected state estimate

    % One-step prediction step: predict the next state and covariance
    x_phat(:, k) = Ad * xhat_k_k + Bd * udev(:, k) + Gd * ddev(:, k); % Predict the next state
    P_k_k1 = Ad * P_k_k * Ad' + Q;       % Predict the next error covariance
    xhat_k_k1 = x_phat(:, k);            % Prepare the next iteration with predicted state
end

end

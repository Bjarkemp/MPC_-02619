function [x_hat, x_phat] = kalman_filter_aug_dynamic_pred(t, xdev, udev, ddev, At, rho, R, Q_aug, Ad_aug, Bd_aug, Gd_aug, Gw_aug, C_aug)
% kalman_filter_aug_dynamic_pred: Implements an augmented dynamic Kalman filter with j-step prediction.
% This function also stores predicted states and covariances for j-step prediction.

% Preallocate output matrices
x_hat = [];                             % Matrix to store corrected state estimates
x_phat = zeros(size(Ad_aug, 1), length(t)); % Matrix to store one-step-ahead state predictions
P_pred = cell(length(t), 100);      % Store predicted error covariance matrices for each step and j-step
x_pred = cell(length(t), 100);      % Store predicted states for each step and j-step

% Cholesky decomposition of measurement noise covariance R
Lr = chol(R, 'lower');                  % Decompose R into lower triangular matrix
v = Lr * (randn(size(xdev)));           % Generate measurement noise with covariance R

% Initialize state estimate and error covariance
xhat_k_k1 = [xdev(:, 1); ddev(:, 1)];   % Initial augmented state estimate
P_k_k1 = 1000 * eye(size(Ad_aug, 1));   % Initial error covariance, large to reflect high uncertainty

% Kalman filter loop over time
for k = 1:length(t)
    % Filtering step: prediction and correction

    % Innovation covariance
    Re_k = C_aug * P_k_k1 * C_aug' + R; % Combine predicted error covariance and measurement noise

    % Kalman gain
    K = P_k_k1 * C_aug' / Re_k;         % Compute optimal Kalman gain

    % Current measurement
    yk = mass_to_height(xdev(:, k), At, rho) + v(:, k); % True measurement with noise

    % Predicted measurement
    yhat_k_k1 = C_aug * xhat_k_k1;      % Predicted measurement using the current estimate

    % Innovation sequence
    ek = yk - yhat_k_k1;                % Difference between actual and predicted measurement

    % Update step: correct the state estimate
    xhat_k_k = xhat_k_k1 + K * ek;      % Update the estimate using the Kalman gain
    P_k_k = P_k_k1 - K * Re_k * K';     % Update the error covariance matrix
    x_hat = [x_hat, xhat_k_k];          % Store the corrected state estimate

    % One-step prediction: predict the next state and error covariance
    x_phat(:, k) = Ad_aug * xhat_k_k + Bd_aug * udev(:, k) + Gd_aug * ddev(:, k); % Predict the next state
    

    % j-step prediction

    P_j = P_k_k1;                       % Start with the current predicted covariance
    for j = 1:100
        % Predict the next state j-steps ahead
        x_phat(:, k+j) = Ad_aug * x_phat(:, k+j-1) + Bd_aug * udev(:, k) + Gd_aug * ddev(:, k);
        % Predict the corresponding covariance
        P_j = Ad_aug * P_j * Ad_aug' + Gw_aug * Q_aug * Gw_aug';
        P_pred{k, j} = P_j;             % Store predicted covariances
    end

    
    xhat_k_k1 = x_phat(:, k);           % Prepare for the next iteration
    % P_k_k1 = Ad_aug * P_k_k * Ad_aug' + Gw_aug * Q_aug * Gw_aug'; % Predict the next error covariance
    P_k_k1 = P_pred{k, end};
end

end

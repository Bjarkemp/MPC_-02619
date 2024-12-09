function [x_hat, x_phat] = kalman_filter_aug_static(t, xdev, udev, ddev, At, rho, R, Q_aug, Ad_aug, Bd_aug, Gd_aug, Gw_aug, C_aug)
x_hat = [];
Lr = chol(R,'lower');                % Cholesky-dekomposition. It just gives me the standard deviation instead of variance.
v = Lr*(randn(size(xdev)));          % Measurement noise. Follows normal distribution with mean=0 and has st.dev of Lr
xhat_k_k1 = [xdev(:,1); ddev(:,1)];  % Start med det initiale tilstandsskøn
P = idare(Ad_aug',C_aug',Gw_aug*Q_aug*Gw_aug',R); % Høj initial kovarians for at tage højde for usikkerhed
for k = 1:length(t)
    % Filtering
    % Innovationskovarians
    Re_k = C_aug * P * C_aug' + R;
    % Kalman-gain
    K = P * C_aug' * inv(Re_k);
    yk = mass_to_height(xdev(:,k),At,rho) + v(:,k);
    yhat_k_k1 = C_aug*xhat_k_k1 + v(:,k);  % Prediktion af måling
    ek = yk - yhat_k_k1;  % Innovationssekvens
    % Update step
    xhat_k_k = xhat_k_k1 + K * ek; 
    % One-step prediction
    x_phat(:,k) = Ad_aug * xhat_k_k + Bd_aug * udev(:,k) + Gd_aug * ddev(:,k);
    % P_k_k = P_k_k1 - K*Re_k*K';
    x_hat = [x_hat xhat_k_k];  % Gem skønnet tilstand
    % Prediktion af fejlkovarians
    % P_k_k1 = Ad_aug * P_k_k * Ad_aug' + Gw_aug*Q_aug*Gw_aug';
    xhat_k_k1 = x_phat(:,k);  % Forbered næste iteration
end


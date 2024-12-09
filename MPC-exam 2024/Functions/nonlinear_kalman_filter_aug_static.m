function [x_hat, x_phat] = nonlinear_kalman_filter_aug_static(t, xdev, udev, ddev,  x0, u, d, At, rho, R, Q_aug, Ad_aug, Gw_aug, C_aug, p)
%--------------------------------------------------------------------------
x_hat = [];
Lr = chol(R,'lower');                % Cholesky-dekomposition. It just gives me the standard deviation instead of variance.
v = Lr*(randn(size(xdev)));          % Measurement noise. Follows normal distribution with mean=0 and has st.dev of Lr
xhat_k_k1 = [xdev(:,1); ddev(:,1)];  % Start med det initiale tilstandsskøn
P = idare(Ad_aug',C_aug',Gw_aug*Q_aug*Gw_aug',R); % Høj initial kovarians for at tage højde for usikkerhed


for k = 1:length(t)-1
    % Filtering
    % Innovationskovarians
    Re_k = C_aug * P * C_aug' + R;
    % Kalman-gain
    K = P * C_aug' * inv(Re_k);
    yk = mass_to_height(xdev(:,k),At,rho) + v(:,k);
    yhat_k_k1 = mass_to_height(xhat_k_k1(1:4),At,rho) + v(:,k);  % Prediktion af måling
    ek = yk - yhat_k_k1;  % Innovationssekvens
    % Update step
    xhat_k_k = xhat_k_k1 + K * ek; 
    % One-step prediction
    x_k_k = xhat_k_k(1:4)+x0;
    [Tk, Xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], x_k_k, ...
        [], u(:,k), d(:,k), p);
    x_phat(:,k) = [Xk(end,:)' - x0; d(:,k) - d(:,1)];

    % P_k_k = P_k_k1 - K*Re_k*K';
    x_hat = [x_hat xhat_k_k];  % Gem skønnet tilstand
    % Prediktion af fejlkovarians
    % P_k_k1 = Ad_aug * P_k_k * Ad_aug' + Gw_aug*Q_aug*Gw_aug';
    xhat_k_k1 = x_phat(:,k);  % Forbered næste iteration
end

k = length(t);

Re_k = C_aug * P * C_aug' + R;
    % Kalman-gain
    K = P * C_aug' * inv(Re_k);
    yk = mass_to_height(xdev(:,k),At,rho) + v(:,k);
    % yhat_k_k1 = C_aug*xhat_k_k1 + v(:,k);  % Prediktion af måling
    yhat_k_k1 = mass_to_height(xhat_k_k1(1:4),At,rho) + v(:,k);  % Prediktion af måling
    ek = yk - yhat_k_k1;  % Innovationssekvens
    % Update step
    xhat_k_k = xhat_k_k1 + K * ek; 

    x_hat = [x_hat xhat_k_k];
%-------------------------------------------------------------------------

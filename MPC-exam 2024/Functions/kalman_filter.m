function [x_hat, x_phat] = kalman_filter(t, xdev, udev, ddev, At, rho, R, Q, Ad, Bd, Gd, C)
                    
Lr = chol(R,'lower');                                   % Cholesky-dekomposition. It just gives me the standard deviation instead of variance.
v = Lr*(randn(size(xdev)));                    % Measurement noise. Follows normal distribution with mean=0 and has st.dev of Lr

xhat_k_k1 = xdev(:,1);  % Start med det initiale tilstandsskøn
P_k_k1 = 10000 * eye(4);  % Høj initial kovarians for at tage højde for usikkerhed

for k = 1:length(t)
    % Filtering

    yhat_k_k1 = C*xhat_k_k1 + v(:,k);  % Prediktion af måling
    ydev(:,k) = C*xdev(:,k) + v(:,k);

    
    ek = ydev(:,k) - yhat_k_k1;  % Innovationssekvens

    % Innovationskovarians
    Re_k = C * P_k_k1 * C' + R;

    % Kalman-gain
    K = P_k_k1 * C' / Re_k;

    % Update step
    xhat_k_k = xhat_k_k1 + K * ek; 
    P_k_k = P_k_k1 - K*Re_k*K';
    % P_k_k = (eye(size(P_k_k1)) - K * C) * P_k_k1;

    x_hat(:,k) = xhat_k_k;  % Gem skønnet tilstand

    % One-step prediction
    x_phat(:,k) = Ad * xhat_k_k + Bd * udev(:,k) + Gd * ddev(:,k);
    xhat_k_k1 = x_phat(:,k);  % Forbered næste iteration

    % Prediktion af fejlkovarians
    P_k_k1 = Ad * P_k_k * Ad' + Q;
end

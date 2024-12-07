function [x_hat, x_phat] = kalman_filter_aug(t, x, u, At, rho, R, Q, Ad, Bd, Gd, C, d)
                    


x_hat = [];
Lr = chol(R,'lower');                                   % Cholesky-dekomposition. It just gives me the standard deviation instead of variance.
v = Lr*(randn(size(xdev)));                    % Measurement noise. Follows normal distribution with mean=0 and has st.dev of Lr

xhat_k_k1 = [xdev(:,1)+3; zeros(2,1)];  % Start med det initiale tilstandsskøn
P_k_k1 = 1 * eye(size(Ad_aug, 1));  % Høj initial kovarians for at tage højde for usikkerhed
for k = 1:length(t)
    % Filtering
    k

    % Innovationskovarians
    Re_k = C_aug * P_k_k1 * C_aug' + R;

    % Kalman-gain
    K = P_k_k1 * C_aug' / Re_k;



    yk = C*xdev(:,k);

    yhat_k_k1 = C_aug*xhat_k_k1 ;  % Prediktion af måling


    ek = yk - yhat_k_k1;  % Innovationssekvens

    % Update step
    xhat_k_k = xhat_k_k1 + K * ek; 

    % One-step prediction
    x_phat(:,k) = Ad_aug * xhat_k_k + Bd_aug * udev(:,k) + Gd_aug * ddev(:,k);


    P_k_k = P_k_k1 - K*Re_k*K';
    % P_k_k = (eye(size(P_k_k1)) - K * C) * P_k_k1;

    x_hat = [x_hat xhat_k_k];  % Gem skønnet tilstand

    % Prediktion af fejlkovarians
    P_k_k1 = Ad_aug * P_k_k * Ad_aug' + Gd_aug*Q_aug*Gd_aug';

    xhat_k_k1 = x_phat(:,k);  % Forbered næste iteration



end


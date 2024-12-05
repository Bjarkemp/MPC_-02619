function [x_hat] = kalman_filter(t, x0, y_sample, At, rho, R, A, B, C)

Lr = chol(R,'lower');

xk_k1 = x0; % Initial state
Pk_k1 = 1000* eye(4); % Initial covariance
x_hat = zeros(4, length(t)); % Pre-allocate estimated states
y_hat = zeros(4, length(t)); % Pre-allocate estimated outputs
x_pred_total =[];
 
for k = 1:length(t)

    % filtering
    yk_k1 = mass_to_height(xk_k1,At,rho);
    ek = y_sample(:,k)-yk_k1;
    Re_k=C*Pk_k1*C'+Lr;
    K = Pk_k1*C'*Re_k^-1;
    xk_k=xk_k1+K*ek; %x_phat(:,k-1) istedet for xk_k1
    Pk_k=Pk_k1-K*Re_k*K';

    x_hat(:,k) = xk_k1; 

    xk_k1 =xk_k;
    Pk_k1 = Pk_k;  

    % One step prediction
    x_phat(:,k) = A*xk_h
end
function [T, X, D, U, x_discrete] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q)

Lq = chol(Q,'lower');                       % Cholesky-dekomposition. It just gives me the standard deviation instead of variance.
w = (Lq*randn(length(x0),length(t)))';              % Measurement noise. Follows normal distribution with mean=0 and has st.dev of Lr

% Initialising compiled masses and time
xk0 = x0;           % [g] Mass at each dicrete time
x_discrete = x0';   % [g] Mass at each dicrete time 
T = [];             % [s] Continuous-time
X = [];             % [g] Mass at all at Continuous-time
D = [];             % [cm3/s] Disturbance for Continuous-time
U = [];             % [cm3/s] Manipulated variables


% Discritized solution of four tank system.
for k=1:length(t)-1
    [Tk, Xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], xk0, ...
        [], u(:,k), d(:,k), p);

    xk0=(Xk(end,:) + w(k,:))';                     % Xk(end:) is Discrete time sample
    x_discrete = [x_discrete; (Xk(end,:) + w(k,:))]; % [g] Mass for Discrete time

    T = [T; Tk];                      
    X = [X; Xk];                      
    D = [D; (d(:,k))'.*ones(length(Xk),2)];
    U = [U; (u(:,k))'.*ones(length(Xk),2)]; 
end
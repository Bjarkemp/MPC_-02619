function [T, X, D, U, x_sample] = discrete_fourtankProcess(x0, t, u, d, p)
% Initialising compiled masses and time
xk0 = x0;           % [g] Mass at time t0
x_sample = x0';     % [g] Sampled mass at time t0
T = [];             % [s] Continuous-time
X = [];             % [g] Mass in all tanks in Continuous-time
D = [];             % [cm3/s] Disturbance in Continuous-time
U = [];             % [cm3/s] Manipulated variables in Continuous-time

% Discritized solution of four tank system.
for k=1:length(t)-1
    [Tk, Xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], xk0, ...
        [], u(:,k), d(:,k), p);

    xk0=Xk(end,:);                        % Xk(end:) is new sample
    x_sample = [x_sample; Xk(end,:)];     % [g] all sampled masses

    T = [T; Tk];                      
    X = [X; Xk];                      
    D = [D; (d(:,k))'.*ones(length(Xk),2)];
    U = [U; (u(:,k))'.*ones(length(Xk),2)]; 
end
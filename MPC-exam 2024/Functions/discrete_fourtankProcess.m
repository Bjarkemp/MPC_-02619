function [T, X, D, U, x_discrete] = discrete_fourtankProcess(x0, t, u, d, p)

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

    xk0=Xk(end,:);                     % Xk(end:) is Discrete time sample
    x_discrete = [x_discrete; Xk(end,:)]; % [g] Mass for Discrete time

    T = [T; Tk];                      
    X = [X; Xk];                      
    D = [D; (d(:,k))'.*ones(length(Xk),2)];
    U = [U; (u(:,k))'.*ones(length(Xk),2)]; 
end
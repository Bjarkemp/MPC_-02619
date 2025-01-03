function [X, U, X_tot, U_tot, T_tot] = PIDcontroller(x0, u0, d, p, N, r, Kc, tau_i, tau_d, t, umin, umax)
% PIDcontroller: Implements a PID controller for the four-tank system.
% The function computes control inputs to regulate tank levels.

% Initialize system state and inputs
x = x0;                             % Initial state (masses in tanks)
z = x0(1:2);                        % Measured tank levels (tanks 1 and 2)
u = u0;                             % Initial input (pump flow rates)
U(:,1) = u0;                        % Store initial input
X(:,1) = x0;                        % Store initial state

% Initialize logging variables
X_tot = [];
U_tot = [];
T_tot = [];

% Initialize controller terms
i = [0; 0];                         % Integral term
e_prev = [0; 0];                    % Previous error for derivative action
dt = (t(end)-t(1))/N;               % Time step size

% Control loop
for k = 1:N
    
    % Simulate system dynamics using the four-tank process
    [T, xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], x, [], u, d(:,k), p);
    
    % Update current state and measured tank levels
    X(:,k+1) = xk(end,:);
    x = xk(end,:);
    z = x(1:2)';

    % ---- PID Controller Calculation --------------------------------------
    e = r(:,k) - z;                            % Compute control error
    i = i + (Kc * dt / tau_i) * e;             % Update integral term
    d_term = (Kc * tau_d / dt) * (e - e_prev); % Update derivative term  
    v = u0 + Kc * e + i + d_term;              % Apply PI control law
    uk = max(umin, min(umax, v));              % Apply input constraints
    % ----------------------------------------------------------------------
    
    % Update control input and store results
    u = uk;
    U(:,k+1) = uk;

    % Update previous error
    e_prev = e;                   

    % Store simulation results
    X_tot = [X_tot; xk];
    U_tot = [U_tot; u'];
    T_tot = [T_tot; T];
    
end
end

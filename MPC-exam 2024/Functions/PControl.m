function [X, U, X_tot, U_tot, T_tot] = Pcontrol(x0, u0, d, p, N, r, Kc, t, umin, umax)
% Pcontrol: Implements a Proportional (P) controller for the four-tank system.

% Initialize state, output, and input variables
x = x0;                          % Current state
z = x0(1:2);                     % Initial measured output (tank levels 1 and 2)
u = u0;                          % Current input
U(:,1) = u0;                     % Store initial input
dk = d(:,1);                     % Initial disturbance
X(:,1) = x0;                     % Store initial state

% Initialize storage for detailed process data
X_tot = [];
U_tot = [];
T_tot = [];

% Main control loop
for k = 1:N
    
    % Simulate the process using the system's dynamics
    [T, xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], x, [], u, d(:,k), p);
    
    % Update and store the last state
    X(:,k+1) = xk(end,:);       
    x = xk(end,:);             
    z = x(1:2)';                      % Update measured output
    
    % ------ P-Controller ---------------------------------------------------
    e = r(:,k) - z;                   % Compute control error
    v = u0 + Kc * e;                  % Apply P control law
    uk = max(umin, min(umax, v));     % Apply input constraints
    % ----------------------------------------------------------------------
    
    % Update input and store data
    u = uk;                      
    U(:,k+1) = uk;               % Store input
    
    % Store detailed trajectories
    X_tot = [X_tot; xk];         % Append state trajectory
    U_tot = [U_tot; u'];         % Append input trajectory
    T_tot = [T_tot; T];          % Append time trajectory
    
end
end

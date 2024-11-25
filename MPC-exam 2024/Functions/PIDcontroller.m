function [X, U, X_tot, U_tot, T_tot] = PIDcontroller(x0, u0, d, p, N, r, Kc, tau_i, tau_d, t, umin, umax)

x = x0;
z = x0(1:2);
u = u0;
U(:,1) = u0;
dk = d(:,1);
X(:,1) = x0;

X_tot = [];
U_tot = [];
T_tot = [];

i = [0; 0]; 
d_term = [0; 0]; 
e_prev = [0; 0]; 
dt = (t(end)-t(1))/N;

for k = 1:N
    
    [T, xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], x, [], u, d(:,k), p);
    X(:,k+1) = xk(end,:);
    x = xk(end,:);
    z = x(1:2)';

    % ------PID-Controller--------------------------------------------------
    e = r-z;
    v = u0 + Kc*e + i + d_term;
    i = i+(Kc*dt/tau_i)*e;
    d_term = (Kc * tau_d / dt) * (e - e_prev);
    uk = max(umin,min(umax,v));
    % ----------------------------------------------------------------------

    u = uk;
    U(:,k+1) = uk;

    e_prev = e;

    X_tot = [X_tot; xk];
    U_tot = [U_tot; u'];
    T_tot = [T_tot; T];
end
function [X, U, X_tot, U_tot, T_tot] = PIcontrol(x0, u0, d, p, N, r, Kc, tau_i, t, umin, umax)

x = x0;
z = x0(1:2);
u = u0;
U(:,1) = u0;
dk = d(:,1);
X(:,1) = x0;

X_tot =[];
U_tot =[];
T_tot =[];

i = [0; 0];
dt = (t(end)-t(1))/N;

for k = 1:N
    [T, xk] = ode15s(@FourTankProcess, [t(k) t(k+1)], x, [], u, d(:,k), p);
    X(:,k+1) = xk(end,:);
    x = xk(end,:);
    z = x(1:2)';

%------PI-Controller-------------------------------------------------------
    e = r-z;
    v = u0 + Kc*e + i;
    i = i+(Kc*dt/tau_i)*e;
    uk = max(umin,min(umax,v));
%--------------------------------------------------------------------------
    u = uk;
    U(:,k+1) = uk;

    X_tot = [X_tot; xk];
    U_tot = [U_tot; u'];
    T_tot = [T_tot; T];
end
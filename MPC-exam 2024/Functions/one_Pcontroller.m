function [X, U] = one_Pcontroller(x0, u0, z0, d, p, N, r, Kc, t, umin, umax)

z = z0;
x = x0;
u = u0;

X(:,1) = x0;
U(:,1) = u0;


for k = 1:N
    % Beregn kontrolinput vha. PI controller
    [u_next] = PControl(r, z, u, Kc, umin, umax);
    u=u_next;
    
    d_step = d(:,k);
    % Løs differentialligningerne for dette tidssteg
    [t_step, x_step] = ode15s(@FourTankProcess, [t(k) t(k+1)], x, [], u, d_step, p);

    % Gem værdier
    X(:,k+1) = x_step(end, :);        % Gem systemtilstandene (højder i tankene)
    U(:,k+1) = u;        % Gem manipulerede variable (flowrater F1 og F2)

    % Opdater tilstanden
    x = x_step(end, :); % Sørg for at x er en kolonnevektor
    z = [x(1); x(2)];
end
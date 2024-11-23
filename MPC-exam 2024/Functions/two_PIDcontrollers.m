function [X, U] = two_PIDcontrollers(x0, u0, z0, d, p, N, i0, r, Kc, taui, t, umin, umax)

z = z0;
x = x0;
u = u0;
i = i0;

X(1,:) = x0;
U(1,:) = u0;

t0 = t(1);
tf = t(end);
dt = (tf-t0)/N;     % [s] interval between each step

for k = 1:N
    for j=1:2
    % Beregn kontrolinput vha. PI controller
        [u_next, i_next] = PIControl(i(j), r(j), z(j), u(j), Kc(j), taui(j), dt, umin(j), umax(j));
        u(j)=u_next;
        i(j)=i_next;
    end
    
    d_step = d(:,k);
    % Løs differentialligningerne for dette tidssteg
    [t_step, x_step] = ode15s(@FourTankProcess, [t(k) t(k+1)], x, [], u, d_step, p);

    % Gem værdier
    X(k+1,:) = x_step(end, :);        % Gem systemtilstandene (højder i tankene)
    U(k+1,:) = u';        % Gem manipulerede variable (flowrater F1 og F2)

    % Opdater tilstanden
    x = x_step(end, :); % Sørg for at x er en kolonnevektor
    z = [x(1) x(2)];
end
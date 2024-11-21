% Define Proportional-Integral Controller
function [u, i] = PIControl(i, r, y, us, Kc, Ti, Ts, umin, umax)
    e = r - y; % Error between setpoint and measured value
    v = us + Kc * e + i; % Control signal with integral action
    i = i + (Kc * Ts / Ti) * e; % Update integral terms
    u = max(umin, min(umax, v)); % Apply control limits
    %u = u(:); % Ensure u is a column vector (2-by-1)
end
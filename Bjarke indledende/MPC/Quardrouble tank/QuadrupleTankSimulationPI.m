
% Main simulation script with controller implementation
t0 = 0;
tf = 1200; % 20 minutes in seconds
Ts = 10; % Sampling time (s)

% Initial parameters
m10 = 0.0; m20 = 0.0; m30 = 0.0; m40 = 0.0;
x0 = [m10; m20; m30; m40];
F1 = 300; F2 = 300;
umin = [0; 0];
umax = [400; 400];
Kc = 0.005; % Controller gain
Ti = 30; % Integral time constant for PI controller

% Setpoints
r = [12000; 10000]; % Desired levels in tanks 1 and 2

% Initial controller states
i = [0; 0]; % Integral terms for PI controller

% Parameters
p = [1.2272; 1.2272; 1.2272; 1.2272; 380.1327; 380.1327; 380.1327; 380.1327; 981; 0.45; 0.40; 1];

% Pre-allocate for performance
num_steps = (tf - t0) / Ts;
X = zeros(num_steps, 4);
T = zeros(num_steps, 1);
U = zeros(num_steps, 2);

% Initial state
x = x0;
t = t0;

for k = 1:num_steps
    % Measurement with noise
    v = 0.5 * randn(4,1); % Noise level
    y = x + v; % Measured levels
    
    % Proportional Controller
    %u = PControl(r, y([1,2]), [F1; F2], Kc, umin, umax);

    % Proportional-Integral Controller
    [u, i] = PIControl(i, r, y([1,2]), [F1; F2], Kc, Ti, Ts, umin, umax);
    
    % Simulate process from t(k) to t(k+1)
    [T_temp, X_temp] = ode15s(@(t,x) QuadrupleTankProcess(t, x, u, p), [t t+Ts], x);
    
    % Update state
    x = X_temp(end, :)';
    T(k) = t;
    X(k, :) = x';
    U(k, :) = u';
    
    % Update time
    t = t + Ts;
end

% Plot results
figure;
subplot(2,1,1);
%plot(T, X./(p(12) * p(5)));
plot(T, X);
xlabel('Time (s)');
%ylabel('Tank liquid height (cm)');
ylabel('Tank liquid mass (m)');
legend('m1', 'm2', 'm3', 'm4');
title('Quadruple Tank System Response with PI Controller');

subplot(2,1,2);
plot(T, U);
xlabel('Time (s)');
ylabel('Control Input (cm^3/s)');
legend('F1', 'F2');
title('Control Inputs over Time');

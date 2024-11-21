% Main simulation script

% Initial parameters
t0 = 0;
tf = 1200; % 20 minutes in seconds
m10 = 0.0; m20 = 0.0; m30 = 0.0; m40 = 0.0;
x0 = [m10; m20; m30; m40];
F1 = 300; F2 = 300;
u = [F1; F2];

% Parameters
p = [1.2272; 1.2272; 1.2272; 1.2272; 380.1327; 380.1327; 380.1327; 380.1327; 981; 0.45; 0.40; 1];

% ODE solver
[T, X] = ode15s(@(t,x) QuadrupleTankProcess(t,x,u,p), [t0 tf], x0);

% Plot results
figure;
plot(T, X);
xlabel('Time (s)');
ylabel('Tank liquid mass (g)');
legend('m1', 'm2', 'm3', 'm4');
title('Quadruple Tank System Response');

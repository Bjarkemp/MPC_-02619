clc; clear all; close all;

% Scalar Standard Brownian Motion = Standard Wiener Process
Ns = 10; % Number of realizations
T = 10; % Final time
N = 1000; % Number of time steps
seed = 100; % Seed for reproducibility
% Realization of Ns Standard Brownian Motions
rng(seed);
dt = T/N;
dW = sqrt(dt)*randn(Ns,N);
W = cumsum(dW,2);

% Time vector
time = linspace(0, T, N);

% Plot the realizations of W(T)
figure;
plot(time, W', 'LineWidth', 1.5); % Plot each realization as a separate line
xlabel('Time (T)');
ylabel('W(T)');
title('Realizations of Standard Brownian Motion');
grid on;
legend(arrayfun(@(x) sprintf('W_%d', x), 1:Ns, 'UniformOutput', false)); % Add legend for each realization
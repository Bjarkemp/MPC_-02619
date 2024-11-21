clc; clear all; close all;

T = 10;
N = 1000;
Ns = 1000;
seed = 100;
[W1,Tw,dW] = ScalarStdWienerProcess(T,N,10,seed);
% Call the Wiener process functions
[W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed);
[Wmean,sW,Wmeanp2sW,Wmeanm2sW] = ScalarSampleMeanStdVar(W);

% Plot realizations of Wiener processes
figure;
plot(Tw, W1, 'linewidth', 2); % Plot all realizations
hold on;

% Plot mean of W
plot(Tw, Wmean, '--r', 'linewidth', 2, 'DisplayName', 'Mean');

% Plot mean + 2*standard deviation
plot(Tw, Wmeanp2sW, '--r', 'linewidth', 2, 'DisplayName', 'Mean + 2*StdDev');

% Plot mean - 2*standard deviation
plot(Tw, Wmeanm2sW, '--r', 'linewidth', 2, 'DisplayName', 'Mean - 2*StdDev');

% Add labels and legend
xlabel('Time (T)');
ylabel('W(t)');
%legend('Realizations', 'Mean', 'Mean + 2 Std Dev', 'Mean - 2 Std Dev');
title('Standard Wiener Process Realizations with Mean and 2*StdDev');
grid off;
hold off;


% x0 = 10;
% lambda = -0.5;
% sigma = 1.0;
% 
% [W,Tw,dW]=ScalarStdWienerProcess(T,N,Ns,seed);
% 
% X = zeros(size(W));
% for i=1:Ns;
% X(i,1) = x0;
% for k=1:N
% dt = Tw(k+1)-Tw(k);
% X(i,k+1) = X(i,k) + lambda*X(i,k)*dt + sigma*dW(i,k);
% end
% end
% 
% plot(Tw, X, 'linewidth', 2); % Plot all realizations
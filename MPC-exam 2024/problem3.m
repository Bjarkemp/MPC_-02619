clc, clear, close all
addpath("Functions");

% Problem 3
% ----------------------------------------------------------
% Parameters
% ----------------------------------------------------------
a1 = 1.2272;                % [cm2] Area of outlet pipe (Tank 1)
a2 = 1.2272;                % [cm2] Area of outlet pipe (Tank 2)
a3 = 1.2272;                % [cm2] Area of outlet pipe (Tank 3)
a4 = 1.2272;                % [cm2] Area of outlet pipe (Tank 4)
A1 = 380.1327;              % [cm2] Cross sectional area (Tank 1)
A2 = 380.1327;              % [cm2] Cross sectional area (Tank 2)
A3 = 380.1327;              % [cm2] Cross sectional area (Tank 3)
A4 = 380.1327;              % [cm2] Cross sectional area (Tank 4)
g = 981;                    % [cm/s2] Gravity
gamma1 = 0.6;               % Flow distribution constant (valve 1)
gamma2 = 0.7;               % Flow distribution constant (valve 2)
rho = 1.0;                  % [g/cm^3] Density of water

% Parameter-vector
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho]; 
at = p(1:4);                % [cm2] Area of outlet pipe
At = p(5:8);                % [cm2] Cross sectional area

% -----------------------------------------------------------
% Simulation scenario
% -----------------------------------------------------------
t0 = 0.0;                   % [s] Initial time
tf= 30*60;                  % [s] End time
dt = 10;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
m10 = 0;                    % [g] Liquid mass in tank 1 at time t0
m20 = 0;                    % [g] Liquid mass in tank 2 at time t0
m30 = 0;                    % [g] Liquid mass in tank 3 at time t0
m40 = 0;                    % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 =100;
F4_0 =150;
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];           % [cm3/s] Disturbance variables at t0
% Brownian motion disturbance
Ns = length(d0); % Number of realizations
seed = 100;
[W,t,dW] = ScalarStdWienerProcess(tf,N,Ns,seed);
sigma = [1^2 0; 0 1^2];        % Covariance for disturbances in F3 and F4
d = d0 + sigma*dW'; 


umin = [0; 0]; % [cm3/s] Manipulated variables
umax = [500; 500]; % [cm3/s] Manipulated variables

h_sp = [25; 30];              % [cm] Setpoint as heigt in tank 1 and tank 2
h_sp_vector = h_sp.*ones(2,length(t));
h_sp_vector(:,80:end) = h_sp_vector(:,80:end)*1.5;
r = h_sp_vector.*(rho * At(1:2));    % [g] Setpoint as mass in tank 1 and tank 2

%% 2 P-controllers with exactly the same parameters

Kc = [0.01]; 

[X1, U1, X_tot1, U_tot1, T_tot1] = Pcontrol(x0, u0, d, p, N, r, Kc, t, umin, umax);
[y1] = sensor_wo_noise(X1', At, rho);

% Plot height of individual tanks
figure(1);
plot(t/60, y1, 'LineWidth', 2);
xlabel('t [min]', 'FontSize', 12);
ylabel('height [cm]', 'FontSize', 12);
legend('T_1', 'T_2', 'T_3', 'T_4', 'Location', 'best');
xlim([0 t(end)/60])
% title('closed loop simulation of liquid height in each tank', 'FontSize', 14);

% Plot height of individual tanks
figure(2);
for i = 1:4
    if i < 3
        subplot(2, 2, i);
        plot(t/60, y1(i,:),'b', 'LineWidth', 2);
        hold on 
        plot(t/60, h_sp_vector(i,:),'k', 'LineWidth', 1);
        hold off
        xlabel('t [min]', 'FontSize', 10);
        ylabel('h [cm]', 'FontSize', 10);
        xlim([0 t(end)/60])
        title(['Tank ', num2str(i)], 'FontSize', 10);
        legend(['h_', num2str(i)], 'Setpoint', 'Location', 'best');
    else
     
    subplot(2, 2, i);
    plot(t/60, y1(i,:),'b', 'LineWidth', 2);
    xlabel('t [min]', 'FontSize', 10);
    ylabel('h [cm]', 'FontSize', 10);
    xlim([0 t(end)/60])
    title(['Tank ', num2str(i)], 'FontSize', 10);
    end
    sgtitle('Liquid Height in each tank', 'FontSize', 12, 'FontWeight', 'bold');
end

% Plot flow rates
figure(3);
for i = 1:2
    subplot(2, 1, i);
    plot(t/60, U1(i,:),'b', 'LineWidth', 2);
    xlabel('t [min]', 'FontSize', 10);
    ylabel('[cm^3/s]', 'FontSize', 10);
    xlim([0 t(end)/60])
    ylim([0 500]);
    title(['F_', num2str(i)], 'FontSize', 10);
end


%% 2 PI-controllers with exactly the same parameters



Kc = [0.01]; 
taui = [2]; % Integral time constant for PI controller

[X2, U2, X_tot2, U_tot2, T_tot2] = PIcontrol(x0, u0, d, p, N, r, Kc, taui, t, umin, umax);
[y2] = sensor_wo_noise(X2', At, rho);

% Plot height of individual tanks
figure(4);
plot(t/60, y2, 'LineWidth', 2);
xlabel('t [min]', 'FontSize', 12);
ylabel('height [cm]', 'FontSize', 12);
legend('T_1', 'T_2', 'T_3', 'T_4', 'Location', 'best');
xlim([0 t(end)/60])
% title('closed loop simulation of liquid height in each tank', 'FontSize', 14);

% Plot height of individual tanks
figure(5);
for i = 1:4
    if i < 3
        subplot(2, 2, i);
        plot(t/60, y2(i,:),'b', 'LineWidth', 2);
        hold on 
        plot(t/60, h_sp_vector(i,:),'k', 'LineWidth', 1);
        hold off
        xlabel('t [min]', 'FontSize', 10);
        ylabel('h [cm]', 'FontSize', 10);
        xlim([0 t(end)/60])
        title(['Tank ', num2str(i)], 'FontSize', 10);
        legend(['h_', num2str(i)], 'Setpoint', 'Location', 'best');
    else
     
    subplot(2, 2, i);
    plot(t/60, y2(i,:),'b', 'LineWidth', 2);
    xlabel('t [min]', 'FontSize', 10);
    ylabel('h [cm]', 'FontSize', 10);
    xlim([0 t(end)/60])
    title(['Tank ', num2str(i)], 'FontSize', 10);
    end
    sgtitle('Liquid Height in each tank', 'FontSize', 12, 'FontWeight', 'bold');
end

% Plot flow rates
figure(6);
for i = 1:2
    subplot(2, 1, i);
    plot(t/60, U2(i,:),'b', 'LineWidth', 2);
    xlabel('t [min]', 'FontSize', 10);
    ylabel('[cm^3/s]', 'FontSize', 10);
    xlim([0 t(end)/60])
    ylim([0 500]);
    title(['F_', num2str(i)], 'FontSize', 10);
end

%% 2 PID-controllers with exactly the same parameters

Kc = [0.01]; 
tau_i = [2]; % Integral time constant for PI controller
tau_d = [80];

[X3, U3, X_tot3 U_tot3, T_tot3] = PIDcontroller(x0, u0, d, p, N, r, Kc, tau_i, tau_d, t, umin, umax);
[y3] = sensor_wo_noise(X3', At, rho);

% Plot height of individual tanks
figure(7);
plot(t/60, y3, 'LineWidth', 2);
xlabel('t [min]', 'FontSize', 12);
ylabel('height [cm]', 'FontSize', 12);
legend('T_1', 'T_2', 'T_3', 'T_4', 'Location', 'best');
xlim([0 t(end)/60])
% title('closed loop simulation of liquid height in each tank', 'FontSize', 14);

% Plot height of individual tanks
figure(8);
for i = 1:4
    if i < 3
        subplot(2, 2, i);
        plot(t/60, y3(i,:),'b', 'LineWidth', 2);
        hold on 
        plot(t/60, h_sp_vector(i,:),'k', 'LineWidth', 1);
        hold off
        xlabel('t [min]', 'FontSize', 10);
        ylabel('h [cm]', 'FontSize', 10);
        xlim([0 t(end)/60])
        title(['Tank ', num2str(i)], 'FontSize', 10);
        legend(['h_', num2str(i)], 'Setpoint', 'Location', 'best');
    else
     
    subplot(2, 2, i);
    plot(t/60, y3(i,:),'b', 'LineWidth', 2);
    xlabel('t [min]', 'FontSize', 10);
    ylabel('h [cm]', 'FontSize', 10);
    xlim([0 t(end)/60])
    title(['Tank ', num2str(i)], 'FontSize', 10);
    end
    sgtitle('Liquid Height in each tank', 'FontSize', 12, 'FontWeight', 'bold');
end

% Plot flow rates
figure(9);
for i = 1:2
    subplot(2, 1, i);
    plot(t/60, U3(i,:),'b', 'LineWidth', 2);
    xlabel('t [min]', 'FontSize', 10);
    ylabel('[cm^3/s]', 'FontSize', 10);
    xlim([0 t(end)/60])
    ylim([0 500]);
    title(['F_', num2str(i)], 'FontSize', 10);
end

%%


% Plot height of individual tanks
figure(10);
for i = 1:4
    if i < 3
        subplot(2, 2, i);
        plot(t/60, y1(i,:),'g', t/60, y2(i,:),'r', t/60, y3(i,:),'b', 'LineWidth', 2);
        hold on 
        plot(t/60, h_sp_vector(i,:),'k', 'LineWidth', 1);
        hold off
        xlabel('t [min]', 'FontSize', 10);
        ylabel('h [cm]', 'FontSize', 10);
        xlim([0 t(end)/60])
        title(['Tank ', num2str(i)], 'FontSize', 10);
        legend(['h_', num2str(i), ' - P-Control'],['h_', num2str(i), ' - PI-Control'],['h_', num2str(i), ' - PID-Control'], 'Setpoint', 'Location', 'best');
    else
     
    subplot(2, 2, i);
        plot(t/60, y1(i,:),'g', t/60, y2(i,:),'r', t/60, y3(i,:),'b', 'LineWidth', 2);
        xlabel('t [min]', 'FontSize', 10);
        ylabel('h [cm]', 'FontSize', 10);
        xlim([0 t(end)/60])
        title(['Tank ', num2str(i)], 'FontSize', 10);
        end
    sgtitle('Liquid Height in each tank', 'FontSize', 12, 'FontWeight', 'bold');
end

% Plot flow rates
figure(9);
for i = 1:2
    subplot(1, 2, i);
    plot(t/60, U1(i,:),'g', t/60, U2(i,:),'r', t/60, U3(i,:),'b', 'LineWidth', 2);
    xlabel('t [min]', 'FontSize', 10);
    ylabel('[cm^3/s]', 'FontSize', 10);
    xlim([0 t(end)/60])
    ylim([0 500]);
    title(['F_', num2str(i)], 'FontSize', 10);
    legend('P-Control','PI-Control','PID-Control', 'Location', 'best');
    sgtitle('Manipulated variables', 'FontSize', 12, 'FontWeight', 'bold');
end
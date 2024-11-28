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
tf= 20*60;                  % [s] End time
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
d = d0.*ones(2, length(t));

umin = [1; 1];
umax = [500; 500];

h_sp = [25; 20];
r = h_sp.*(rho * At(1:2)); 

%% 2 P-controllers with exactly the same parameters

Kc = [0.05]; 

[X1, U1, X_tot1, U_tot1, T_tot1] = Pcontrol(x0, u0, d, p, N, r, Kc, t, umin, umax);
[y1] = sensor_wo_noise(X1', At, rho);

% Plot height of individual tanks
figure(1);
plot(t/60, y1, 'LineWidth', 2);
xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('\textbf{[cm]}', 'FontSize', 12, 'Interpreter', 'latex');
legend('T_1', 'T_2', 'T_3', 'T_4', 'Location', 'best');
xlim([0 t(end)/60])
% title('closed loop simulation of liquid height in each tank', 'FontSize', 14);

% Plot height of individual tanks
figure(2);
for i = 1:4
    subplot(2, 2, i);
    plot(t/60, y1(i,:), 'LineWidth', 2);
    xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('\textbf{h [cm]}', 'FontSize', 12, 'Interpreter', 'latex');
    xlim([0 t(end)/60])
    % ylim([0 50]);
    title(['Height in tank ', num2str(i)], 'FontSize', 10);
end

% Plot flow rates
figure(3);
for i = 1:2
    subplot(1, 2, i);
    plot(t/60, U1(i,:), 'LineWidth', 2);
    xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('Flow [cm^3/s]', 'FontSize', 12);
    xlim([0 t(end)/60])
    ylim([0 500]);
    title(['F', num2str(i)], 'FontSize', 10);
end


%% 2 PI-controllers with exactly the same parameters

Kc = [0.05]; 
taui = [12]; % Integral time constant for PI controller

[X2, U2, X_tot2, U_tot2, T_tot2] = PIcontrol(x0, u0, d, p, N, r, Kc, taui, t, umin, umax);
[y2] = sensor_wo_noise(X2', At, rho);

% Plot height of individual tanks
figure(4);
plot(t/60, y2, 'LineWidth', 2);
xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('\textbf{[cm]}', 'FontSize', 12, 'Interpreter', 'latex');
legend('T_1', 'T_2', 'T_3', 'T_4', 'Location', 'best');
xlim([0 t(end)/60])
% title('closed loop simulation of liquid height in each tank', 'FontSize', 14);

% Plot height of individual tanks
figure(5);
for i = 1:4
    subplot(2, 2, i);
    plot(t/60, y2(i,:), 'LineWidth', 2);
    xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('\textbf{h [cm]}', 'FontSize', 12, 'Interpreter', 'latex');
    xlim([0 t(end)/60])
    ylim([0 50]);
    title(['Height in tank ', num2str(i)], 'FontSize', 10);
end

% Plot flow rates
figure(6);
for i = 1:2
    subplot(1, 2, i);
    plot(t/60, U2(i,:), 'LineWidth', 2);
    xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('Flow [cm^3/s]', 'FontSize', 12);
    xlim([0 t(end)/60])
    ylim([0 500]);
    title(['F', num2str(i)], 'FontSize', 10);
end

%% 2 PID-controllers with exactly the same parameters

Kc = [0.05]; 
tau_i = [12]; % Integral time constant for PI controller
tau_d = [20];

[X3, U3, X_tot3 U_tot3, T_tot3] = PIDcontroller(x0, u0, d, p, N, r, Kc, tau_i, tau_d, t, umin, umax);
[y3] = sensor_wo_noise(X3', At, rho);

% Plot height of individual tanks
figure(4);
plot(t/60, y3, 'LineWidth', 2);
xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('\textbf{[cm]}', 'FontSize', 12, 'Interpreter', 'latex');
legend('T_1', 'T_2', 'T_3', 'T_4', 'Location', 'best');
xlim([0 t(end)/60])
% title('closed loop simulation of liquid height in each tank', 'FontSize', 14);

% Plot height of individual tanks
figure(5);
for i = 1:4
    subplot(2, 2, i);
    plot(t/60, y3(i,:), 'LineWidth', 2);
    xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('\textbf{h [cm]}', 'FontSize', 12, 'Interpreter', 'latex');
    xlim([0 t(end)/60])
    ylim([0 50]);
    title(['Height in tank ', num2str(i)], 'FontSize', 10);
end

% Plot flow rates
figure(6);
for i = 1:2
    subplot(1, 2, i);
    plot(t/60, U3(i,:), 'LineWidth', 2);
    xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('Flow [cm^3/s]', 'FontSize', 12);
    xlim([0 t(end)/60])
    ylim([0 500]);
    title(['F', num2str(i)], 'FontSize', 10);
end
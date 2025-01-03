clc, clear, close all
addpath("Functions");

% Problem 4
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
dt = 1;                    % [s] interval between each step
N = tf/dt;                  % Number of steps 
t = t0:dt:tf;               % [s] time-vector
m10 = 17612.0123864868;                    % [g] Liquid mass in tank 1 at time t0
m20 = 29640.6694933624;                    % [g] Liquid mass in tank 2 at time t0
m30 = 4644.21948249842;                    % [g] Liquid mass in tank 3 at time t0
m40 = 9378.49308238599;                    % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;                 % [cm3/s] Flow rate from pump 1
F2_0 = 300;                 % [cm3/s] Flow rate from pump 2
F3_0 = 100;
F4_0 = 150;
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [F3_0; F4_0;];           % [cm3/s] Disturbance variables at t0
d = d0.*ones(2, length(t));
u = u0.*ones(2, length(t));
[y0] = sensor_wo_noise(x0', At, rho);

%%  -------------------- 4.1 10% step change ------------------------------
% Stepchange in u1

u_stepchange = find(t==5*60); % Time that step change occur
stepchange = 1 + 0.1;         % Step change size
                              % Adding step change to u1 only
u(1,u_stepchange:end) = u0(1)*stepchange; 

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y_10_u1] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

% Deviation variables
ydev = y_10_u1-y0;                  
udev = u-u0;

% plots(t,udev,ydev')
% sgtitle('Samlet Titel for Plottet');
% hold off
% 
% % Steady state gain for T1 and T2
% K1 = ydev(1:2,end) / udev(1,end);
% 
% % finding time constants (assuming first order)
% tau1_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
% tau1_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
% 
% tau1 = [tau1_1; tau1_2];
% 
% %genrate transfer functions for 1st order system
% sys11 = tfest(udev',ydev(1:2,:)',1)
% % sys = tfest(udev',ydev',1,'Ts',dt);
% % sys_continuous = d2c(sys)
% 
% %generate transfer fucktions for 2nd order system
% sys21 = tfest(udev',ydev(1:2,:)',2)
% % sys = tfest(udev',ydev',1,'Ts',dt);
% % sys_continuous = d2c(sys)

%% Step change in u2
clc
% Reset the manipulated variables to the start value
u = u0.*ones(2, length(t));

% Adding step change to u2 only
u(2,u_stepchange:end) = u0(2)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y_10_u2] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

% Deviation variables
ydev = y_10_u2-y0;
udev = u-u0;

% plots(t,udev,ydev')
% sgtitle('Samlet Titel for Plottet');
% hold off
% 
% % Steady state gain for T1 and T2
% K2 = ydev(1:2,end) / udev(2,end);
% 
% % finding time constants (assuming second order)
% tau2_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
% tau2_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
% 
% tau2 = [tau2_1; tau2_2];
% 
% %genrate transfer functions for 1st order system
% sys12 = tfest(udev',ydev(1:2,:)',1)
% % sys = tfest(udev',ydev',1,'Ts',dt);
% % sys_continuous = d2c(sys)
% %generate transfer fucktions for 2nd order system
% sys22 = tfest(udev',ydev(1:2,:)',2)
% % sys = tfest(udev',ydev',1,'Ts',dt);
% % sys_continuous = d2c(sys)

%%  -------------------- 4.1 25% step change ------------------------------
close all

% Reset the manipulated variables
u = u0.*ones(2, length(t));

u_stepchange = find(t==5*60);
stepchange = 1 + 0.25;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y_25_u1] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y_25_u1-y0;
udev = u-u0;

% plots(t,udev,ydev')
% sgtitle('Samlet Titel for Plottet');
% hold off
% 
% K1 = ydev(:,end) / udev(1,end);
% % finding time constants
% tau1_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
% tau1_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
% tau1_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end))));
% tau1_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end)))) - t(u_stepchange);
% 
% tau1 = [tau1_1; tau1_2; tau1_3; tau1_4];


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y_25_u2] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y_25_u2-y0;
udev = u-u0;
% 
% plots(t,udev,ydev')
% sgtitle('Samlet Titel for Plottet');
% hold off
% 
% K2 = ydev(:,end) / udev(2,end);
% % finding time constants
% tau2_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
% tau2_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
% tau2_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end)))) - t(u_stepchange);
% tau2_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end))));
% 
% tau2 = [tau2_1; tau2_2; tau2_3; tau2_4];

%%  -------------------- 4.1 50% step change ------------------------------

% Reset the manipulated variables
u = u0.*ones(2, length(t));

u_stepchange = find(t==5*60);
stepchange = 1+0.5;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y_50_u1] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y_50_u1-y0;
udev = u-u0;

% plots(t,udev,ydev')
% sgtitle('Samlet Titel for Plottet');
% hold off

K1 = ydev(:,end) / udev(1,end);
% finding time constants
tau1_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau1_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
tau1_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end))));
tau1_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end)))) - t(u_stepchange);

tau1 = [tau1_1; tau1_2; tau1_3; tau1_4];


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y_50_u2] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y_50_u2-y0;
udev = u-u0;

% plots(t,udev,ydev')
% sgtitle('Samlet Titel for Plottet');
% hold off

K2 = ydev(:,end) / udev(2,end);
% finding time constants
tau2_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau2_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
tau2_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end)))) - t(u_stepchange);
tau2_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end))));

tau2 = [tau2_1; tau2_2; tau2_3; tau2_4];

% Convert time vector to minutes
t_minutes = t / 60;

% Simulated responses for Tanks 1 to 4 under 10%, 25%, and 50% step changes in u1 and u2
y1_10_u1 = y_10_u1(1, :); % Tank 1 response, 10% step in u1
y1_25_u1 = y_25_u1(1, :); % Tank 1 response, 25% step in u1
y1_50_u1 = y_50_u1(1, :); % Tank 1 response, 50% step in u1

y2_10_u1 = y_10_u1(2, :); % Tank 2 response, 10% step in u1
y2_25_u1 = y_25_u1(2, :); % Tank 2 response, 25% step in u1
y2_50_u1 = y_50_u1(2, :); % Tank 2 response, 50% step in u1

y3_10_u1 = y_10_u1(3, :); % Tank 3 response, 10% step in u1
y3_25_u1 = y_25_u1(3, :); % Tank 3 response, 25% step in u1
y3_50_u1 = y_50_u1(3, :); % Tank 3 response, 50% step in u1

y4_10_u1 = y_10_u1(4, :); % Tank 4 response, 10% step in u1
y4_25_u1 = y_25_u1(4, :); % Tank 4 response, 25% step in u1
y4_50_u1 = y_50_u1(4, :); % Tank 4 response, 50% step in u1

% Repeat for u2 step changes
y1_10_u2 = y_10_u2(1, :); % Tank 1 response, 10% step in u2
y1_25_u2 = y_25_u2(1, :); % Tank 1 response, 25% step in u2
y1_50_u2 = y_50_u2(1, :); % Tank 1 response, 50% step in u2

y2_10_u2 = y_10_u2(2, :); % Tank 2 response, 10% step in u2
y2_25_u2 = y_25_u2(2, :); % Tank 2 response, 25% step in u2
y2_50_u2 = y_50_u2(2, :); % Tank 2 response, 50% step in u2

y3_10_u2 = y_10_u2(3, :); % Tank 3 response, 10% step in u2
y3_25_u2 = y_25_u2(3, :); % Tank 3 response, 25% step in u2
y3_50_u2 = y_50_u2(3, :); % Tank 3 response, 50% step in u2

y4_10_u2 = y_10_u2(4, :); % Tank 4 response, 10% step in u2
y4_25_u2 = y_25_u2(4, :); % Tank 4 response, 25% step in u2
y4_50_u2 = y_50_u2(4, :); % Tank 4 response, 50% step in u2

% Plot for step changes in u1
figure;
subplot(2, 2, 1);
hold on;
plot(t_minutes, y1_10_u1, '-r', 'LineWidth', 1.5);
plot(t_minutes, y1_25_u1, '-b', 'LineWidth', 1.5);
plot(t_minutes, y1_50_u1, '-k', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 1');
legend({'10%', '25%', '50%'}, 'Location', 'northwest');

subplot(2, 2, 2);
hold on;
plot(t_minutes, y2_10_u1, '-r', 'LineWidth', 1.5);
plot(t_minutes, y2_25_u1, '-b', 'LineWidth', 1.5);
plot(t_minutes, y2_50_u1, '-k', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 2');
legend({'10%', '25%', '50%'}, 'Location', 'northwest');

subplot(2, 2, 3);
hold on;
plot(t_minutes, y3_10_u1, '-r', 'LineWidth', 1.5);
plot(t_minutes, y3_25_u1, '-b', 'LineWidth', 1.5);
plot(t_minutes, y3_50_u1, '-k', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 3');
legend({'10%', '25%', '50%'}, 'Location', 'northwest');

subplot(2, 2, 4);
hold on;
plot(t_minutes, y4_10_u1, '-r', 'LineWidth', 1.5);
plot(t_minutes, y4_25_u1, '-b', 'LineWidth', 1.5);
plot(t_minutes, y4_50_u1, '-k', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 4');
legend({'10%', '25%', '50%'}, 'Location', 'northwest');

sgtitle('Step Responses for u_1 Changes');

%% save figure
saveas(gcf,fullfile('C:\Users\bjark\OneDrive\Skrivebord\MPC_-02619\MPC-exam 2024\Plots','problem4u1.png'),'png')

% Plot for step changes in u2
figure;
subplot(2, 2, 1);
hold on;
plot(t_minutes, y1_10_u2, '-r', 'LineWidth', 1.5);
plot(t_minutes, y1_25_u2, '-b', 'LineWidth', 1.5);
plot(t_minutes, y1_50_u2, '-k', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 1');
legend({'10%', '25%', '50%'}, 'Location', 'northwest');

subplot(2, 2, 2);
hold on;
plot(t_minutes, y2_10_u2, '-r', 'LineWidth', 1.5);
plot(t_minutes, y2_25_u2, '-b', 'LineWidth', 1.5);
plot(t_minutes, y2_50_u2, '-k', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 2');
legend({'10%', '25%', '50%'}, 'Location', 'northwest');

subplot(2, 2, 3);
hold on;
plot(t_minutes, y3_10_u2, '-r', 'LineWidth', 1.5);
plot(t_minutes, y3_25_u2, '-b', 'LineWidth', 1.5);
plot(t_minutes, y3_50_u2, '-k', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 3');
legend({'10%', '25%', '50%'}, 'Location', 'northwest');

subplot(2, 2, 4);
hold on;
plot(t_minutes, y4_10_u2, '-r', 'LineWidth', 1.5);
plot(t_minutes, y4_25_u2, '-b', 'LineWidth', 1.5);
plot(t_minutes, y4_50_u2, '-k', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 4');
legend({'10%', '25%', '50%'}, 'Location', 'northwest');

sgtitle('Step Responses for u_2 Changes');

%% save figure
saveas(gcf,fullfile('C:\Users\bjark\OneDrive\Skrivebord\MPC_-02619\MPC-exam 2024\Plots','problem4u2.png'),'png')

%%  -------------------- 4.2 Low noise 25% step change --------------------

close all

u = u0.*ones(2, length(t));

R = [0.4^2 0 0 0; 0 0.5^2 0 0; 0 0 0.05^2 0; 0 0 0 0.1^2];     % Covariance for measurement noise
Q = [40^2 0 0 0; 0 50^2 0 0; 0 0 5^2 0; 0 0 0 10^2];     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.25;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y_low] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

%%  -------------------- 4.2 medium noise 25% step change --------------------

close all

u = u0.*ones(2, length(t));

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*2;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*2;     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.25;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y_med] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

%%  -------------------- 4.2 high noise 25% step change --------------------

close all

u = u0.*ones(2, length(t));

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.25;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y_high] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

% Simulated responses for Tanks 1 to 4 under 10% step changes in u1 with
% noise
y1_low = y_low(1, :); % Tank 1 response, 10% step in u1
y1_med = y_med(1, :); % Tank 1 response, 25% step in u1
y1_high = y_high(1, :); % Tank 1 response, 50% step in u1

y2_low = y_low(2, :); % Tank 2 response, 10% step in u1
y2_med = y_med(2, :); % Tank 2 response, 25% step in u1
y2_high = y_high(2, :); % Tank 2 response, 50% step in u1

y3_low = y_low(3, :); % Tank 3 response, 10% step in u1
y3_med = y_med(3, :); % Tank 3 response, 25% step in u1
y3_high = y_high(3, :); % Tank 3 response, 50% step in u1

y4_low = y_low(4, :); % Tank 4 response, 10% step in u1
y4_med = y_med(4, :); % Tank 4 response, 25% step in u1
y4_high = y_high(4, :); % Tank 4 response, 50% step in u1

% Plot for step changes in u2
figure;
subplot(2, 2, 1);
hold on;
plot(t_minutes, y1_high, '-k', 'LineWidth', 1.5);
plot(t_minutes, y1_med, '-b', 'LineWidth', 1.5);
plot(t_minutes, y1_low, '-r', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 1');
legend({'low', 'medium', 'high'}, 'Location', 'southeast');

subplot(2, 2, 2);
hold on;
plot(t_minutes, y2_high, '-k', 'LineWidth', 1.5);
plot(t_minutes, y2_med, '-b', 'LineWidth', 1.5);
plot(t_minutes, y2_low, '-r', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 2');
legend({'low', 'medium', 'high'}, 'Location', 'southeast');

subplot(2, 2, 3);
hold on;
plot(t_minutes, y3_high, '-k', 'LineWidth', 1.5);
plot(t_minutes, y3_med, '-b', 'LineWidth', 1.5);
plot(t_minutes, y3_low, '-r', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 3');
legend({'low', 'medium', 'high'}, 'Location', 'southeast');

subplot(2, 2, 4);
hold on;
plot(t_minutes, y4_high, '-k', 'LineWidth', 1.5);
plot(t_minutes, y4_med, '-b', 'LineWidth', 1.5);
plot(t_minutes, y4_low, '-r', 'LineWidth', 1.5);
hold off;
xlabel('Time [min]');
ylabel('Height [cm]');
title('Tank 4');
legend({'low', 'medium', 'high'}, 'Location', 'southeast');

sgtitle('Step Responses for 10% u_1 Change with different noise');

%% save figure
saveas(gcf,fullfile('C:\Users\bjark\OneDrive\Skrivebord\MPC_-02619\MPC-exam 2024\Plots','problem4noise.png'),'png')

%% transfer

clear all; close all; clc;

% Define transfer functions
% u1 -> y1: 0.0016 / (s + 0.0101)
sys_u1_y1 = tf(0.0016, [1, 0.0101]);

% u1 -> y2: (0.0618e-4 * s + 0.1479e-4) / (s^2 + 0.022 * s + 0.0001)
sys_u1_y2 = tf([0.0618e-4, 0.1479e-4], [1, 0.022, 0.0001]);

% u2 -> y1: (0.0544e-4 * s + 0.1571e-4) / (s^2 + 0.0305 * s + 0.0002)
sys_u2_y1 = tf([0.0544e-4, 0.1571e-4], [1, 0.0305, 0.0002]);

% u2 -> y2: 0.0018 / (s + 0.0079)
sys_u2_y2 = tf(0.0018, [1, 0.0079]);

% Time vector for step response
t = 0:1:20*60; % Time in seconds

% Compute step responses
[y_u1_y1, t_u1_y1] = step(sys_u1_y1, t);
[y_u1_y2, t_u1_y2] = step(sys_u1_y2, t);
[y_u2_y1, t_u2_y1] = step(sys_u2_y1, t);
[y_u2_y2, t_u2_y2] = step(sys_u2_y2, t);

% Plot step responses
figure;

% Plot u1 -> y1
subplot(2, 2, 1);
plot(t_u1_y1, y_u1_y1, 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Response');
title('u_1 \rightarrow y_1');
grid on;
xlim([0 20*60])
% Plot u1 -> y2
subplot(2, 2, 2);
plot(t_u1_y2, y_u1_y2, 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Response');
title('u_1 \rightarrow y_2');
grid on;
xlim([0 20*60])
% Plot u2 -> y1
subplot(2, 2, 3);
plot(t_u2_y1, y_u2_y1, 'g', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Response');
title('u_2 \rightarrow y_1');
grid on;
xlim([0 20*60])
% Plot u2 -> y2
subplot(2, 2, 4);
plot(t_u2_y2, y_u2_y2, 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Response');
title('u_2 \rightarrow y_2');
grid on;
xlim([0 20*60])

% Add a global title for the figure
sgtitle('Step Responses of Identified Transfer Functions');

saveas(gcf,fullfile('C:\Users\bjark\OneDrive\Skrivebord\MPC_-02619\MPC-exam 2024\Plots','problem4tf.png'),'png')

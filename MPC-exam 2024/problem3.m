clc, clear, close all
addpath("Functions");
%% Problem 3
% ----------------------------------------------------------
% Parameters
% ----------------------------------------------------------
a1 = 1.2272;   % [cm2] Area of outlet pipe (Tank 1)
a2 = 1.2272;   % [cm2] Area of outlet pipe (Tank 2)
a3 = 1.2272;   % [cm2] Area of outlet pipe (Tank 3)
a4 = 1.2272;   % [cm2] Area of outlet pipe (Tank 4)
A1 = 380.1327; % [cm2] Cross sectional area (Tank 1)
A2 = 380.1327; % [cm2] Cross sectional area (Tank 2)
A3 = 380.1327; % [cm2] Cross sectional area (Tank 3)
A4 = 380.1327; % [cm2] Cross sectional area (Tank 4)
g = 981;       % [cm/s2] Gravity
gamma1 = 0.6;  % Flow distribution constant (valve 1)
gamma2 = 0.7;  % Flow distribution constant (valve 2)
rho = 1.0;     % [g/cm^3] Density of water

% Parameter-vector
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho]; 
at = p(1:4);   % [cm2] Area of outlet pipe
At = p(5:8);   % [cm2] Cross sectional area

% -----------------------------------------------------------
% Simulation scenario
% -----------------------------------------------------------
t0 = 0.0;           % [s] Initial time
tf= 20*60;          % [s] End time
N = 200;             % Number of steps 
dt = (tf-t0)/N;     % [s] interval between each step
t = t0:dt:tf;       % [s] time-vector
m10 = 0;            % [g] Liquid mass in tank 1 at time t0
m20 = 0;            % [g] Liquid mass in tank 2 at time t0
m30 = 0;            % [g] Liquid mass in tank 3 at time t0
m40 = 0;            % [g] Liquid mass in tank 4 at time t0
F1_0 = 300;         % [cm3/s] Flow rate from pump 1
F2_0 = 300;         % [cm3/s] Flow rate from pump 2
x0 = [m10; m20; m30; m40];    % [g] Start values 
u0 = [F1_0; F2_0];            % [cm3/s] Manipulated variables 
d0 = [0; 0;];                 % [cm3/s] Disturbance variables at t0
F_0 = [u0',d0'];              % [cm3/s] MV & DV at t0
% u = u0.*ones(2, length(t));
d = d0.*ones(2, length(t));



% Notice that the valve position determines if the system is reachable or
% not. It is not certain that the desired set points can be reached if the
% valves has a certain position that does not allow it.


% minimum & maximum values of manipulated variables
umin = [1 1];
umax = [500 500];

Kc = [0.002 0.002]; % Controller gain
taui = [2000 2000]; % Integral time constant for PI controller

h_sp = [30 20]; % Height Set point

% Converts height set point to mass reference
r = h_sp.*(rho*p(5:6)'); 

% Initial controller states
i = [0 0]; % Integral terms for PI controller. The accumilated error is 0 to t0


%Initialising mass, time and manipulated variables
X = zeros(N, 4);
T = zeros(N, 1);
U = zeros(N, 2);

x = x0;
z = [x0(1) x0(2)];
u = [u0(1) u0(2)];

%%
%% PI control loop
for k = 1:N

    for j=1:2
    % Beregn kontrolinput vha. PI controller
        [u_temp, i_temp] = PIControl(i(j), r(j), z(j), u(j), Kc(j), taui(j), dt, umin(j), umax(j));
        u(j)=u_temp;
        i(j)=i_temp;
    end
    % Løs differentialligningerne for dette tidssteg
    [T_step, X_step] = ode15s(@FourTankProcess, [t0:dt:tf], x, [], u, d, p);

    % Gem værdier
    T(k) = t(k);          % Gem tiden
    X(k,:) = x;        % Gem systemtilstandene (højder i tankene)
    X(k,:) = X(k,:) + X(k,:).*(0.05*randn(1,length(X(k,:))));
    U(k,:) = u';        % Gem manipulerede variable (flowrater F1 og F2)
    U(k,:) = U(k,:) + U(k,:).*(0.05*randn(1,length(U(k,:))));

    % Opdater tilstanden
    x = X_step(end, :); % Sørg for at x er en kolonnevektor
    z = [x(1) x(2)];
end

%%
% Helper function to compute height
convert_mass_to_height = @(X, rho, p) X ./ (rho * p(5:8)');
% Convert mass to height
    h = convert_mass_to_height(X, rho, p);



%% Plot
figure(1);
plot(T/60, h, 'LineWidth', 2);
grid on;
xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('\textbf{[cm]}', 'FontSize', 12, 'Interpreter', 'latex');
legend('T_1', 'T_2', 'T_3', 'T_4', 'Location', 'best');
% title('closed loop simulation of liquid height in each tank', 'FontSize', 14);


% Plot height of individual tanks
figure(2);
for i = 1:4
    subplot(2, 2, i);
    plot(T/60, h(:, i), 'LineWidth', 2);
    xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('\textbf{h [cm]}', 'FontSize', 12, 'Interpreter', 'latex');
    % xlim([0 T_total(end)]);
    title(['Height in tank ', num2str(i)], 'FontSize', 10);
    grid on;
end

% Plot flow rates
figure(3);
for i = 1:2
    subplot(1, 2, i);
    plot(T/60, U(:, i), 'LineWidth', 2);
    xlabel('\textbf{t [min]}', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('Flow [cm^3/s]', 'FontSize', 12);
    % xlim([0 T_total(end)]);
    ylim([0 500]);
    title(['F', num2str(i)], 'FontSize', 10);
    grid on;
end
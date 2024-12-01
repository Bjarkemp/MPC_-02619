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
F3_0 =100;
F4_0 =150;
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

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

% Deviation variables
ydev = y-y0;                  
udev = u-u0;

plots(t,udev,ydev')
sgtitle('Samlet Titel for Plottet');
hold off

% Steady state gain for T1 and T2
K1 = ydev(1:2,end) / udev(1,end);

% finding time constants (assuming first order)
tau1_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau1_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);

tau1 = [tau1_1; tau1_2];

%genrate transfer functions for 1st order system
sys11 = tfest(udev',ydev(1:2,:)',1)
% sys = tfest(udev',ydev',1,'Ts',dt);
% sys_continuous = d2c(sys)

%generate transfer fucktions for 2nd order system
sys21 = tfest(udev',ydev(1:2,:)',2)
% sys = tfest(udev',ydev',1,'Ts',dt);
% sys_continuous = d2c(sys)

%% Step change in u2
clc
% Reset the manipulated variables to the start value
u = u0.*ones(2, length(t));

% Adding step change to u2 only
u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

% Deviation variables
ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')
sgtitle('Samlet Titel for Plottet');
hold off

% Steady state gain for T1 and T2
K2 = ydev(1:2,end) / udev(2,end);

% finding time constants (assuming second order)
tau2_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau2_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);

tau2 = [tau2_1; tau2_2];

%genrate transfer functions for 1st order system
sys12 = tfest(udev',ydev(1:2,:)',1)
% sys = tfest(udev',ydev',1,'Ts',dt);
% sys_continuous = d2c(sys)
%generate transfer fucktions for 2nd order system
sys22 = tfest(udev',ydev(1:2,:)',2)
% sys = tfest(udev',ydev',1,'Ts',dt);
% sys_continuous = d2c(sys)

%%  -------------------- 4.1 25% step change ------------------------------
close all

% Reset the manipulated variables
u = u0.*ones(2, length(t));

u_stepchange = find(t==5*60);
stepchange = 1 + 0.25;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')
sgtitle('Samlet Titel for Plottet');
hold off

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

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')
sgtitle('Samlet Titel for Plottet');
hold off

K2 = ydev(:,end) / udev(2,end);
% finding time constants
tau2_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau2_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
tau2_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end)))) - t(u_stepchange);
tau2_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end))));

tau2 = [tau2_1; tau2_2; tau2_3; tau2_4];

%%  -------------------- 4.1 50% step change ------------------------------

% Reset the manipulated variables
u = u0.*ones(2, length(t));

u_stepchange = find(t==5*60);
stepchange = 1+0.5;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')
sgtitle('Samlet Titel for Plottet');
hold off

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

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')
sgtitle('Samlet Titel for Plottet');
hold off

K2 = ydev(:,end) / udev(2,end);
% finding time constants
tau2_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau2_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
tau2_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end)))) - t(u_stepchange);
tau2_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end))));

tau2 = [tau2_1; tau2_2; tau2_3; tau2_4];

%%  -------------------- 4.2 Low noise 10% step change --------------------
close all

u = u0.*ones(2, length(t));

R = [0.4^2 0 0 0; 0 0.5^2 0 0; 0 0 0.05^2 0; 0 0 0 0.1^2];     % Covariance for measurement noise
Q = [40^2 0 0 0; 0 50^2 0 0; 0 0 5^2 0; 0 0 0 10^2];     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.1;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')
sgtitle('Samlet Titel for Plottet');
hold off


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')
sgtitle('Samlet Titel for Plottet');
hold off

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

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')

%%  -------------------- 4.2 Low noise 50% step change --------------------

close all

u = u0.*ones(2, length(t));

R = [0.4^2 0 0 0; 0 0.5^2 0 0; 0 0 0.05^2 0; 0 0 0 0.1^2];     % Covariance for measurement noise
Q = [40^2 0 0 0; 0 50^2 0 0; 0 0 5^2 0; 0 0 0 10^2];     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.5;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


%%  -------------------- 4.2 medium noise 10% step change --------------------
close all

u = u0.*ones(2, length(t));

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*2;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*2;     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.1;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


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

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')

%%  -------------------- 4.2 medium noise 50% step change --------------------

close all

u = u0.*ones(2, length(t));

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*2;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*2;     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.5;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')

%%  -------------------- 4.2 high noise 10% step change --------------------
close all

u = u0.*ones(2, length(t));

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.1;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


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

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')

%%  -------------------- 4.2 high noise 50% step change --------------------

close all

u = u0.*ones(2, length(t));

R = [(0.4)^2 0 0 0; 0 (0.5)^2 0 0; 0 0 (0.05)^2 0; 0 0 0 (0.1)^2]*4;     % Covariance for measurement noise
Q = [(40)^2 0 0 0; 0 (50)^2 0 0; 0 0 (5)^2 0; 0 0 0 (10)^2]*4;     % Covariance for process noise

u_stepchange = find(t==5*60);
stepchange = 1 + 0.5;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*stepchange;


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')

%%  ------------------------------ 4.3 ------------------------------------

% Her beregnes hvad masse ved steady state er i med det givne u0, d0 og p. 
% Jeg regnede det fordi jeg gerne ville bruge det som startværdi, så det
% ikke behøvede at stabiliserer sig før vi laver step change.

xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);    % Løser differentialligningssystemet
% (QuadrupleTankProcess) og solver hvad x er når hældningen er 0. Dvs. at
% den beregner hvad masserne er når der opnås steady state i tankene.

%%
%genrate transfer functions for 1st order system
sys1 = tfest(udev',ydev(1:2,:)',2)
% sys = tfest(udev',ydev',1,'Ts',dt);
% sys_continuous = d2c(sys)
<<<<<<< Updated upstream
%generate transfer fucktions for 2nd order system
sys2 = tfest(udev',ydev(1:2,:)',2)
% sys = tfest(udev',ydev',1,'Ts',dt);
% sys_continuous = d2c(sys)

%% state space model
% Extract the cells and turn to matrix for change in u1
denominators11 = sys11.denominator; 
aU11= cell2mat(denominators11(1:2, 1));
numerators11 = sys11.Numerator;
bY11 = cell2mat(numerators11(1:2,1));
denominators21 = sys21.denominator; 
aU21= cell2mat(denominators21(1:2, 1));
numerators21 = sys21.Numerator;
bY21 = cell2mat(numerators21(1:2,1));

% Extract the cells and turn to matrix for change in u2
denominators12 = sys12.denominator; 
aU12= cell2mat(denominators12(1:2, 1));
numerators12 = sys12.Numerator;
bY12 = cell2mat(numerators12(1:2,1));
denominators22 = sys22.denominator; 
aU22= cell2mat(denominators22(1:2, 1));
numerators22 = sys22.Numerator;
bY22 = cell2mat(numerators22(1:2,1));

%u1y1
b11 = bY11(1,:);
a11 = aU11(1,:);
[A11,B11,C11,D11] = tf2ss(b11,a11);

%u1y2
b21 = bY21(2,:);
a21 = aU21(2,:);
[A12,B12,C12,D12] = tf2ss(b21,a21);

%u2y1
b22 = bY22(1,:);
a22 = aU22(1,:);
[A21,B21,C21,D21] = tf2ss(b22,a22);

%y2u2
b12 = bY11(1,:);
a12 = aU11(1,:);
[A22,B22,C22,D22] = tf2ss(b12,a12);
=======




%% 4.6 
>>>>>>> Stashed changes

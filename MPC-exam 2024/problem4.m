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
dt = 10;                    % [s] interval between each step
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

u_stepchange = find(t==5*60);
stepchange = 1 + 0.1;
u(1,u_stepchange:end) = u0(1)*stepchange;

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')

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

K2 = ydev(:,end) / udev(2,end);
% finding time constants
tau2_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau2_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
tau2_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end)))) - t(u_stepchange);
tau2_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end))));

tau2 = [tau2_1; tau2_2; tau2_3; tau2_4];

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

R = [(2*0.4)^2 0 0 0; 0 (2*0.5)^2 0 0; 0 0 (2*0.05)^2 0; 0 0 0 (2*0.1)^2];     % Covariance for measurement noise
Q = [(2*40)^2 0 0 0; 0 (2*50)^2 0 0; 0 0 (2*5)^2 0; 0 0 0 (2*10)^2];     % Covariance for process noise

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

%%  -------------------- 4.2 medium noise 50% step change --------------------

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

%%  -------------------- high noise ------------------------------------

u = u0.*ones(2, length(t));

R = [(4*0.4)^2 0 0 0; 0 (4*0.5)^2 0 0; 0 0 (4*0.05)^2 0; 0 0 0 (4*0.1)^2];     % Covariance for measurement noise
Q = [(4*40)^2 0 0 0; 0 (4*50)^2 0 0; 0 0 (4*5)^2 0; 0 0 0 (4*10)^2];     % Covariance for process noise

u_stepchange1 = find(t==40*60);
u(1,u_stepchange1:end) = u0(1)*(1-0.1);

u_stepchange2 = find(t==80*60);
u(1,u_stepchange2:end) = u(1,end)*(1+0.25);

u_stepchange3 = find(t==120*60);
u(1,u_stepchange3:end) = u(1,end)*(1-0.5);


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y] = sensor_plus_noise(x, At, rho, R);
[z] = output(x, At, rho);
plots(t,u,y')


u = u0.*ones(2, length(t));

u_stepchange1 = find(t==40*60);
u(2,u_stepchange1:end) = u0(2)*(1-0.1);

u_stepchange2 = find(t==80*60);
u(2,u_stepchange2:end) = u(2,end)*(1+0.25);

u_stepchange3 = find(t==120*60);
u(2,u_stepchange3:end) = u(2,end)*(1-0.5);


% Solve ODE for this step
[T2, X2, D2, U2, x2] = discrete_fourtankProcess_plus_noise(x0, t, u, d, p, Q);

[y2] = sensor_plus_noise(x2, At, rho, R);
[z2] = output(x2, At, rho);
plots(t,u,y2')

%%  ------------------------------ 4.3 ------------------------------------

% Her beregnes hvad masse ved steady state er i med det givne u0, d0 og p. 
% Jeg regnede det fordi jeg gerne ville bruge det som startværdi, så det
% ikke behøvede at stabiliserer sig før vi laver step change.

xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);    % Løser differentialligningssystemet
% (QuadrupleTankProcess) og solver hvad x er når hældningen er 0. Dvs. at
% den beregner hvad masserne er når der opnås steady state i tankene.


%% Steady state gains

% For first step cahnge
K1_1 = (y1s2-y1s1)/(u1s2-u1s1);
K2_1 = (y2s2-y2s1)/(u2s2-u2s1);
% For second step cahnge
K1_2 = (y1s3-y1s2)/(u1s3-u1s2);
K2_2 = (y2s3-y2s2)/(u2s3-u2s2);
% For third step cahnge
K1_3 = (y1s4-y1s3)/(u1s4-u1s3);
K2_3 = (y2s4-y2s3)/(u2s4-u2s3);

% Compiling the steady state gains for each step change
K1_all = [K1_1 K1_2 K1_3];
K2_all = [K2_1 K2_2 K2_3];

% Min forventning var faktisk at ligemeget hvilken stepchange, så burde
% systemet have samme steady state gains, men de er faktisk en smule
% forskellige og det synes jeg er mærkeligt

%%  -------------------- 4.1 ---------------------------------------------

[y0] = sensor_wo_noise(x0', At, rho);

u_stepchange = find(t==5*60);
u(1,u_stepchange:end) = u0(1)*(1+0.1);

% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')

K1 = ydev(:,end) / udev(1,end);
% finding time constants
tau1_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau1_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
tau1_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end))));
tau1_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end)))) - t(u_stepchange);

tau1 = [tau1_1; tau1_2; tau1_3; tau1_4];


% Reset the manipulated variables
u = u0.*ones(2, length(t));

u(2,u_stepchange:end) = u0(2)*(1+0.1);


% Solve ODE for this step
[T, X, D, U, x] = discrete_fourtankProcess(x0, t, u, d, p);

[y] = sensor_wo_noise(x, At, rho);
[z] = output(x, At, rho);

ydev = y-y0;
udev = u-u0;

plots(t,udev,ydev')

K2 = ydev(:,end) / udev(2,end);
% finding time constants
tau2_1 = min(t(find(ydev(1,:)>=0.632*ydev(1,end)))) - t(u_stepchange);
tau2_2 = min(t(find(ydev(2,:)>=0.632*ydev(2,end)))) - t(u_stepchange);
tau2_3 = min(t(find(ydev(3,:)>=0.632*ydev(3,end)))) - t(u_stepchange);
tau2_4 = min(t(find(ydev(4,:)>=0.632*ydev(4,end))));

tau2 = [tau2_1; tau2_2; tau2_3; tau2_4];



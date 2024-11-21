function xdot = QuadrupleTankProcess(t,x,u,p)
% -----------------------------------------------------------
% Parameters
% -----------------------------------------------------------
a1 = 1.2272; %[cm2] Area of outlet pipe 1
a2 = 1.2272; %[cm2] Area of outlet pipe 2
a3 = 1.2272; %[cm2] Area of outlet pipe 3
a4 = 1.2272; %[cm2] Area of outlet pipe 4
A1 = 380.1327; %[cm2] Cross sectional area of tank 1
A2 = 380.1327; %[cm2] Cross sectional area of tank 2
A3 = 380.1327; %[cm2] Cross sectional area of tank 3
A4 = 380.1327; %[cm2] Cross sectional area of tank 4
g = 981; %[cm/s2] The acceleration of gravity
gamma1 = 0.45; % Flow distribution constant. Valve 1
gamma2 = 0.40; % Flow distribution constant. Valve 2
rho = 1; %[g/cm2] Density of liquid
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2];

% ------------------------------------------------------------
% Simulation scenario
% ------------------------------------------------------------
t0 = 0.0; % [s] Initial time
tf = 20*60; % [s] Final time
h10 = 0.0; % [cm] Liquid level in tank 1 at time t0
h20 = 0.0; % [cm] Liquid level in tank 2 at time t0
h30 = 0.0; % [cm] Liquid level in tank 3 at time t0
h40 = 0.0; % [cm] Liquid level in tank 4 at time t0
F1 = 300; % [cm3/s] Flow rate from pump 1
F2 = 300; % [cm3/s] Flow rate from pump 2
x0 = [h10; h20; h30; h40];
u = [F1; F2];

% ------------------------------------------------------------
% Equations
% ------------------------------------------------------------

m1 = rho * A1 * h1;
m2 = rho * A2 * h2;
m3 = rho * A3 * h3;
m4 = rho * A4 * h4;

q1in = gamma1 * F1;
q2in = gamma2 * F2;
q3in = (1-gamma2) * F2;
q4in = (1-gamma1) * F1;

q1 = a1 * sqrt(2*g*h1);
q2 = a2 * sqrt(2*g*h2);
q3 = a3 * sqrt(2*g*h3);
q4 = a4 * sqrt(2*g*h4);

%diff

dm1dt = rho * q1in + rho * q3 - rho * q1;
dm2dt = rho * q2in + rho * q4 - rho * q2;
dm3dt = rho * q3in - rho * q3;
dm4dt = rho * q4in - rho * q4;

% Measurement and Output Function
y(:,k) = g(x(:,k)) + v(:,k);
z(:,k) = h(x(:,k));
% Controller
u(:,k) = Controller(r(:,k),y(:,k),d(:,k));
% Simulate process from t(k) to t(k+1)
[T,X]=ode15s(@f, [t(k) t(k+1)], x(:,k), ...
odeOptions, ...
u(:,k), d(:,k), w(:,k) );
x(:,k+1) = X(end,:)';

[T,X]=ode15s(@QuadrupleTankProcess,...
[t0 tf], x0, ODEoptions, u, p)
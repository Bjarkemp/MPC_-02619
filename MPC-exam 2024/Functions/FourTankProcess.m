function dxdt = FourTankProcess(t,x,u,d,p)
a = p(1:4);                        % Pipe cross sectional areas [cm^2]
A = p(5:8);                        % Tank cross sectional areas [cm^2]
g = p(9);                          % Acceleration of gravity [cm/s^2]
gamma = p(10:11);                  % [-] Valve positions 
rho = p(12);                       % [g/cm^3] Density of water 

F3 = d(1);                         % [cm3/s] Disturbance to tank 3
F4 = d(2);                         % [cm3/s] Disturbance to tank 4
m = x;                             % [g] Liquid mass in each tank
  
h = m./(rho*A);                    % [cm] Liquid level in each tank 
q = a.*sqrt(2*g*h);                % [cm3/s] Outflow from each tank 

qin=zeros(4,1);                    % Initializing inflows
qin(1,1) = gamma(1)*u(1);          % [cm3/s] Inflow to tank 1
qin(2,1) = gamma(2)*u(2);          % [cm3/s] Inflow to tank 2
qin(3,1) = (1-gamma(2))*u(2);      % [cm3/s] Inflow to tank 3
qin(4,1) = (1-gamma(1))*u(1);      % [cm3/s] Inflow to tank 4

% Initialize time derivatives of mass in each tank
dxdt=zeros(4,1);                   
dxdt(1,1)=rho*(qin(1)+q(3)-q(1));  % [g/s] Rate of change of mass Tank 1
dxdt(2,1)=rho*(qin(2)+q(4)-q(2));  % [g/s] Rate of change of mass Tank 2
dxdt(3,1)=rho*(qin(3)-q(3)+F3);    % [g/s] Rate of change of mass Tank 3
dxdt(4,1)=rho*(qin(4)-q(4)+F4);    % [g/s] Rate of change of mass Tank 4
end
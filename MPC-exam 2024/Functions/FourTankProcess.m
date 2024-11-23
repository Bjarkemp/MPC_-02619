function dxdt = FourTankProcess(t,x,u, d, p)
%This function implements a differential equation model for the 4-tank system.

 m=x;               % Mass of liquid in each tank [g]
 F=u;               % Flow rates in pumps [cm3/s]
 a=p(1:4);        % Pipe cross sectional areas [cm2]
 A=p(5:8);        % Tank cross sectional areas [cm2]
 g=p(9,1);          % Valve positions [-]
 gamma=p(10:11,1);  % Acceleration of gravity [cm/s2]
 rho=p(12,1);       % Density of water [g/cm3]
  
 % Compute height
 h=m./(rho*A);      % Liquid level in each tank [cm]

 % Compute outflows
 q=a.*sqrt(2*g*h);  % Outflow from each tank [cm3/s]

 %compute inflows
 qin=zeros(4,1);
 qin(1,1)=gamma(1)*F(1);
 qin(2,1)=gamma(2)*F(2);
 qin(3,1)=(1-gamma(2))*F(2);
 qin(4,1)=(1-gamma(1))*F(1);

% Compute dxdt
dxdt=zeros(4,1);
dxdt(1,1)=rho*(qin(1)+q(3)-q(1));
dxdt(2,1)=rho*(qin(2)+q(4)-q(2));
dxdt(3,1)=rho*(qin(3)-q(3)+d(1));
dxdt(4,1)=rho*(qin(4)-q(4)+d(2));
 end
clc, clear, close all
% addpath("Functions");

% Problem 5
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
dt = 20;                    % [s] interval between each step
N = 10;                     % Number of steps 
t = t0:dt:tf;               % [s] time-vector
Ph = 5;                     % Prediction horizon
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

%%  -------------------- 8.1 MPC function ------------------------------

%linearization
% Steady State
xs = fsolve(@FourTankSystemWrap,x0,[],u0,d0,p);
ys = sensor_wo_noise(xs,at,rho);
zs = sensor_wo_noise(xs,at,rho);

%Stochastic Brownian
R = [1^2 0 0 0; 0 1^2 0 0; 0 0 0.5^2 0; 0 0 0 0.5^2];     % Covariance for measurement noise
sigma = [2^2 0; 0 2^2]; 
[A, B, C, G, Gw] = linearized_models(xs, At, at, rho, gamma1, gamma2, g, ...
                                                                                'brownian', d0, sigma, R);
%-------------------------------------------------------------------------
% Construct the matrix for exponential calculation
M = [-A, G * G';
     zeros(size(A')), A'];

% Compute the matrix exponential
Phi = expm(M * dt);

% Extract submatrices
Phi_11 = Phi(1:size(A,1), 1:size(A,1));
Phi_12 = Phi(1:size(A,1), size(A,1)+1:end);
Phi_22 = Phi(size(A,1)+1:end, size(A,1)+1:end);

% Compute process noise covariance matrix Q
Q2 = Phi_22' * Phi_12;

%ZOH Discretization of Linear System
%Stochastic Brownian
[Ad,Bd,Gd]=c2dzoh(A,B,G,dt);
% sys = ss(A_k,[B_k,E_k],C_c,[D_c,zeros(size(E_c,2))]);
D = zeros(2,4);
sys = ss(Ad,[Bd,Gd],C(1:2,:),D);

C = sys.C;
Q_hat = Q2;
[sys_aug,Q_hat_aug] = sys_aug_func(sys,Q_hat);

Ts = 50;
Q = 100;% Tuning parameter
S = eye(2)*0.01;% Tuning parameter
t_R = 10*60; % The reference is constant for 10 minutes
t_f = t_R*5; % New time in step reference
N = t_R/Ts; % Should be at the same length at time horizon of the dynamics

MPC_sys = MPCDesign(sys,Q,S,N);

% noise
sigma_R = 0.5; % COVARIANCE
Rd = 5;
v_k = mvnrnd(zeros(t_f/Ts,size(C,1)),eye(2)*sigma_R)';
d_k = mvnrnd(zeros(t_f/Ts,size(Gd,2)),eye(2)*Rd)';

x = zeros(4,t_f/Ts+1);
u_vec = zeros(2,t_f/Ts+1);
R1 = [0 ; 0];
R2 = [10 ; 10];
R3 = [20 ;-10];

R=zeros(2*(N+t_f/Ts),1);
R(1:2*(t_f/Ts/3),1)=kron(ones(t_f/Ts/3,1),R1);
R(2*(t_f/Ts/3)+1:2*(2*t_f/Ts/3),1)=kron(ones(t_f/Ts/3,1),R2);
R(2*(2*t_f/Ts/3)+1:2*(N+t_f/Ts),1)=kron(ones(N+t_f/Ts/3,1),R3);

x_s = x0';
inputs = {x,x_s,u_vec,R,v_k,d_k};

cons = []; % Constraints

type = 1; % Unconstrained MPC

[y,y_nl,u_vec,u_nl] = MPC_Sim(sys,sys_aug,MPC_sys,[],Q_hat,Q_hat_aug,sigma_R,t_f,t_R,Ts,inputs,N,cons,[],[],p,type);

% Data used to plot
t_plot = (1:Ts:t_f)/60;
R_plot = [ones(t_R/Ts,1)*R1(1)+h_s(1) , ones(t_R/Ts,1)*R1(2)+h_s(2) ;
               ones(t_R/Ts,1)*R2(1)+h_s(1) , ones(t_R/Ts,1)*R2(2)+h_s(2) ;
               ones(t_R/Ts,1)*R3(1)+h_s(1) , ones(t_R/Ts,1)*R3(2)+h_s(2)];

figure
subplot (2,1,1); hold on
plot(t_plot,y(1,:)+h_s(1),'r');plot(t_plot,y_nl(1,:),'--m'); plot(t_plot,R_plot(:,1),'--k'); plot(t_plot,y(2,:)+h_s(2),'b');plot(t_plot,y_nl(2,:),'--g'); plot(t_plot,R_plot(:,2),'-.k')
hold off; grid on; legend('h_1 - lin','h_1 - Nonlinear system','r_1','h_2 - lin','h_2 - Nonlinear system','r_2','Location','best'); ylabel('Height $[cm]$','Interpreter','latex')
subplot (2,1,2); hold on
plot(t_plot,u_vec(1,:)+u(1),'r');plot(t_plot,u_nl(1,:)+u(1),'--m'); plot(t_plot,u_vec(2,:)+u(2),'b');plot(t_plot,u_nl(2,:)+u(2),'--g');
hold off; grid on; legend('u_1','u_1 - Nonlinear system','u_2','u_2 - Nonlinear system','Location','best'); xlabel('Time [min]'), ylabel('Flow $[cm^3/s]$')

function MPC_sys = MPCDesign(sys,Q,S,N)
    [phi,phi_w,Gamma,Gamma_d] = MPC_Constants(sys,N);
    % Determine Qz weight
    Q_z = eye(size(Gamma))*Q;
    % Generate Hs
    Hs = zeros(N*2);
    % Create arrow for matrix dimensions
    pil_Hs_row = [1:size(S,1)]; pil_Hs_col = pil_Hs_row;
    for i = 1:N
        while (pil_Hs_row(end) < N*2)
            Hs(pil_Hs_row,pil_Hs_col) = S*2; % Set S in diagonal
            Hs(pil_Hs_row,pil_Hs_col+size(S,1)) = -S;
            Hs(pil_Hs_row+size(S,1),pil_Hs_col) = -S;
            % Update arrows
            pil_Hs_row = pil_Hs_row + size(S,1); pil_Hs_col = pil_Hs_row; 
        end % while
        if i == N
            Hs(pil_Hs_row,pil_Hs_col) = S;
        end % if i == N
        % Update arrows
        pil_Hs_row = [1:size(S,1)] + size(S,1)*i; pil_Hs_col = pil_Hs_row;
    end % i
    
    % Determine H
    H = (transpose(Gamma) * Q_z * Gamma + Hs);
    H = (H + H') / 2; % Ensure symmetry

    
    % Now determine elements for g
    M_x0 = transpose(Gamma)*Q_z*phi;
    M_r = -transpose(Gamma)*Q_z;
    M_d  = transpose(Gamma)*Q_z*Gamma_d;
    M_u1 = zeros(size(M_d,1),size(S,1)); M_u1(1:size(S,1),1:size(S,1)) = -S;

    % Now we determine the bounds on the input according to
    % "Lecture_07C_MPC" slide 18
    Lambda = eye(N);
    Lambda = kron(eye(2),Lambda);
    for i=1:N*2
        Lambda(i+2:i+3,i:i+1) = -eye(2);
    end
    Lambda = Lambda(1:end-3,1:end-1);
    I_0 = zeros(N*size(S,1),size(S,2));
    I_0(1:size(S,1),1:size(S,2)) = eye(size(S,1));
    
    MPC_sys = {phi,phi_w,Gamma,Gamma_d,Q_z,M_x0,M_r,M_d,M_u1,Hs,H,Lambda,I_0};
end % function

function [phi,phi_w,Gamma,Gamma_d] = MPC_Constants(sys,N)
    % Determine matrice
    A = sys.A; B = sys.B(:,1:2); E = sys.B(:,3:4);C = sys.C;
    
    % First create arrows to ensure matrix dimensions
    pil_phi = [1:size(C,1)];
    Nul_mat = zeros(size(C,1),2);
    pil_Gamma_row = [1:size(C*A*B,1)];
    pil_Gamma_col = [1:size(C*A*B,2)];
    pil_Gamma_d_row = [1:size(C*A*E,1)];
    pil_Gamma_d_col = [1:size(C*A*E,2)];
    
    for i = 1:N
        % Determine phi for N iterations.
        phi(pil_phi,:) = C*A^i;
        phi_w(pil_phi,:) = C*A^i*E;
        pil_phi = pil_phi+size(C,1);

        % Determine Gamma for N iterations
        while(pil_Gamma_row(end)< N*size(Nul_mat,2)+1)
            Gamma(pil_Gamma_row,pil_Gamma_col) = [C*A^(i-1)*B];
            pil_Gamma_row = pil_Gamma_row + size(Nul_mat,2);
            pil_Gamma_col = pil_Gamma_col + size(Nul_mat,1);
        end% while
        % Shift row and columns
        pil_Gamma_row = [1:size(Nul_mat,1)] + size(Nul_mat,1)*i;
        pil_Gamma_col = [1:size(Nul_mat,2)];
        
        % Determine Gamma_d for N iterations
        while(pil_Gamma_d_row(end)< N*size(Nul_mat,2)+1)
            Gamma_d(pil_Gamma_d_row,pil_Gamma_d_col) = [C*A^(i-1)*E];
            pil_Gamma_d_row = pil_Gamma_d_row + size(Nul_mat,2);
            pil_Gamma_d_col = pil_Gamma_d_col + size(Nul_mat,1);
        end % while
        % Shift row and columns
        pil_Gamma_d_row = [1:size(Nul_mat,1)] + size(Nul_mat,1)*i;
        pil_Gamma_d_col = [1:size(Nul_mat,2)];
    end % i
end

function [y,y_nl,u,u_nl] = MPC_Sim(sys,sys_aug,MPC_sys,MPC_sys_soft,Q_hat,Q_hat_aug,sigma_R,t_f,t_R,Ts,inputs,N,cons,cons_InOut,cons_Econ,p,type)
    x = cell2mat(inputs(1,1));  x_s = cell2mat(inputs(1,2)); u = cell2mat(inputs(1,3)); 
    R = cell2mat(inputs(1,4)); v_k = cell2mat(inputs(1,5)); d_k = cell2mat(inputs(1,6));
    h_s = ( x_s(:,1:2)/(p(5)*p(12)) )';
    
    % System
    A = sys.A; B = sys.B(:,1:2); E = sys.B(:,3:4);C = sys.C;
    A_aug = sys_aug.A; B_aug = sys_aug.B(:,1:2); E_aug = sys_aug.B(:,3:4);C_aug = sys_aug.C;
    
    % Initialize
    x_hat(:,1) = [x(:,1) ; zeros(size(E,2),1)];
    x_nl = x; x_nl(:,1) = x_s';
    x_hat_nl(:,1) = zeros(6,1);
    u_nl = u;
    
    % Kalman 
    R_noise = eye(size(C,1),size(C,1)) * sqrt(sigma_R);
    P = idare(A_aug',C_aug',Q_hat_aug,R_noise);
    Re = C_aug*P*C_aug'+R_noise;
    L = P*C_aug'*inv(Re);
        
    t_nl = 0;
for i = 1:t_f/Ts
    if i == 1 % ONLY FIRST ITERATION
        if type == 1 % UNCONSTRAINED
            u_mpc = Uncon_MPC(MPC_sys,x_hat(1:4,i),R(1:2*N,1),u(:,i));
            u_mpc_nl = Uncon_MPC(MPC_sys,x_hat_nl(1:4,i),R(1:2*N,1),u_nl(:,i));
        elseif type ==2 % INPUT CONSTRAINED
            u_mpc = U_con_MPC(MPC_sys,x_hat(1:4,i),R(1:2*N,1),u(:,i),cons);
            u_mpc_nl = U_con_MPC(MPC_sys,x_hat_nl(1:4,i),R(1:2*N,1),u_nl(:,i),cons);
        elseif type ==3 % INPUT OUTPUT CONSTRAINED
            u_mpc =  InOut_MPC(MPC_sys,MPC_sys_soft,x_hat(1:4,i),R(1:2*N,1),u(:,i),cons_InOut,N);
            u_mpc_nl =  InOut_MPC(MPC_sys,MPC_sys_soft,x_hat_nl(1:4,i),R(1:2*N,1),u_nl(:,i),cons_InOut,N);
        elseif type ==4 % Economic MPC
            u_mpc = Econ_MPC(MPC_sys,x_hat(1:4,i),R(1:2*N,1),u(:,i),cons,cons_Econ);
            u_mpc_nl = Econ_MPC(MPC_sys,x_hat_nl(1:4,i),R(1:2*N,1),u_nl(:,i),cons,cons_Econ);
        end
    else % end i == 1 
        if type == 1% UNCONSTRAINED
            u_mpc = Uncon_MPC(MPC_sys,x_hat(1:4,i),R(2*i+1:2*N+2*i,1),u(:,i-1));
            u_mpc_nl = Uncon_MPC(MPC_sys,x_hat_nl(1:4,i),R(2*i+1:2*N+2*i,1),u_nl(:,i-1));
        elseif type == 2% INPUT CONSTRAINED
            u_mpc = U_con_MPC(MPC_sys,x_hat(1:4,i),R(2*i+1:2*N+2*i,1),u(:,i-1),cons);
            u_mpc_nl = U_con_MPC(MPC_sys,x_hat_nl(1:4,i),R(2*i+1:2*N+2*i,1),u_nl(:,i-1),cons);
        elseif type == 3 % INPUT OUTPUT CONSTRAINED
            u_mpc =  InOut_MPC(MPC_sys,MPC_sys_soft,x_hat(1:4,i),R(2*i+1:2*N+2*i,1),u(:,i-1),cons_InOut,N);
            u_mpc_nl =  InOut_MPC(MPC_sys,MPC_sys_soft,x_hat_nl(1:4,i),R(2*i+1:2*N+2*i,1),u_nl(:,i-1),cons_InOut,N);
        elseif type ==4 % Economic MPC
            u_mpc = Econ_MPC(MPC_sys,x_hat(1:4,i),R(2*i+1:2*N+2*i,1),u(:,i-1),cons,cons_Econ);
            u_mpc_nl = Econ_MPC(MPC_sys,x_hat_nl(1:4,i),R(2*i+1:2*N+2*i,1),u_nl(:,i-1),cons,cons_Econ);
        end
    end % i
    u(:,i) = u_mpc(1:2,1); % Can only use from 1 to number of inputs
    u_nl(:,i) = u_mpc_nl(1:2,1); % Can only use from 1 to number of inputs
    
    % Non linear sim
    F_nl = [u_nl(:,i)'+[300 300] d_k(:,i)'+[250 250]];
    [t,x_temp] = ode15s(@FourTankSystem,[t_nl(i) t_nl(i)+Ts],x_nl(1:4,i),[],F_nl,p);
    t_nl(i+1) = t(end);
    x_nl(:,i+1) = x_temp(end,:)';
    y_nl(:,i) = x_nl(1:2,i)/(p(12)*p(5)) + v_k(:,i);
    
    y_hat_nl(:,i) = C_aug*x_hat_nl(:,i);
    e_nl(:,i) = (y_nl(:,i) - h_s) - y_hat_nl(:,i);
    x_hat_nl(:,i+1) = A_aug*x_hat_nl(:,i) + B_aug*u_nl(:,i)  +  L*e_nl(:,i);
    
    % Linear sim
    x(:,i+1) = A*x(:,i) + B*u(:,i) + E*d_k(:,i);
    y(:,i) = C*x(:,i) + v_k(:,i);   
    
    y_hat(:,i) = C_aug*x_hat(:,i);
    e(:,i) = y(:,i) - y_hat(:,i);        
    x_hat(:,i+1) = A_aug*x_hat(:,i) + B_aug*u(:,i)  +  L*e(:,i);       
end % i      
u = u(:,1:end-1);
u_nl = u_nl(:,1:end-1); 
end % function

function [sys_aug,Q_hat_aug] = sys_aug_func(sys,Q_hat)
    A = sys.A; B = sys.B(:,1:2); E = sys.B(:,3:4);C = sys.C; D = sys.D(1:2,1:2);

    A_aug = [A E ; zeros(size(E,2),size(A,1)) eye(size(E,2))]; % Aug with E insted of B!!
    B_aug = [B ; zeros(size(B,2))];
    C_aug = [C , zeros(size(C,1))];
    E_aug = [E ; zeros(size(E,2))];
        
    Q_hat_aug = [Q_hat , zeros(size(Q_hat,1),size(E,2)) ; zeros(size(E,2),size(Q_hat,1)) , eye(size(E,2))];

    sys_aug = ss(A_aug,[B_aug,E_aug],C_aug,[D,zeros(size(E,2))]);
end

function u_mpc = Uncon_MPC(MPC_sys,x,R,u_prev)
    % Extract data for design
    M_x0 = cell2mat(MPC_sys(1,6)); M_r = cell2mat(MPC_sys(1,7)); M_u1 = cell2mat(MPC_sys(1,9));
    Hs = cell2mat(MPC_sys(1,10)); H = cell2mat(MPC_sys(1,11));
    g = M_x0*x + M_r*R + M_u1*u_prev;
    [u_mpc] = qpsolver(H+Hs,g,[],[],[],[],[],[]);
end % function

function [x_new] = qpsolver(H,g,low_x,up_x,mat,lb,ub,x_init)
    A_quad = [mat ; -mat];
    bound = [ub ; -lb];
    options=optimoptions("quadprog","Display","none");
    [x_new info] = quadprog(H,g,A_quad,bound,[],[],low_x,up_x,x_init,options);
end

function xdot = FourTankSystem(t,x,u,p)
% FOURTANKSYSTEM Model dx/dt = f(t,x,u,p) for 4-tank System
%
% This function implements a differential equation model for the
% 4-tank system.
%
% Syntax: xdot = FourTankSystem(t,x,u,p)
% Unpack states, MVs, and parameters
m = x;              % Mass of liquid in each tank [g]
F = u;              % Flow rates in pumps [cm3/s]
a = p(1:4,1);       % Pipe cross sectional areas [cm2]
A = p(5:8,1);       % Tank cross sectional areas [cm2]
gamma = p(9:10,1);  % Valve positions [-]
g = p(11,1);        % Acceleration of gravity [cm/s2]
rho = p(12,1);      % Density of water [g/cm3]

% Inflows
qin = zeros(4,1);
qin(1,1) = gamma(1)*F(1);     % Inflow from valve 1 to tank 1 [cm3/s]
qin(2,1) = gamma(2)*F(2);     % Inflow from valve 2 to tank 2 [cm3/s]
qin(3,1) = (1-gamma(2))*F(2); % Inflow from valve 2 to tank 3 [cm3/s]
qin(4,1) = (1-gamma(1))*F(1); % Inflow from valve 1 to tank 4 [cm3/s]

% Outflows
h = m./(rho*A);             % Liquid level in each tank [cm]
qout = a.*sqrt(2*g*h);      % Outflow from each tank [cm3/s]

% Differential equations
xdot = zeros(4,1);
xdot(1,1) = rho*(qin(1,1)+qout(3,1)-qout(1,1)); % Mass balance Tank 1
xdot(2,1) = rho*(qin(2,1)+qout(4,1)-qout(2,1)); % Mass balance Tank 2
xdot(3,1) = rho*(qin(3,1)-qout(3,1)+F(3));      % Mass balance Tank 3
xdot(4,1) = rho*(qin(4,1)-qout(4,1)+F(4));      % Mass balance Tank 4
end
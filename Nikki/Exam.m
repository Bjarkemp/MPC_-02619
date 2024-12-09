clc; clear all; close all;

t_f = 20*60;

%% Problem 3.1 Step Deterministic
F_ss = [u' d'];
Per = [1.1 ; 1.25 ; 1.5];

% Generate 3-dimensional matrix with steadystate inputs, and steps on each
% input
F_in = []; F_in = repmat(F_ss,[3,1]);
for i = 1:2
    F_in = cat(3,F_in,repmat(F_ss,[3,1]));
    for j = 1:3
        F_per = F_ss*Per(j);
        F_in(j,i,i+1) = F_per(i);
    end % j
end % i 

% Run initla steady state simulation
[T,X] = ode15s(@FourTankSystem,[0 t_f],x0,[],F_in(1,:,1),p);
x_s = X(end,:);
t_0 = T(end)';
t_end = t_f+t_0;

n_subplot = 4; % Number of subplots
StateNames = ['Tank 3'; 'Tank 4';'Tank 1';'Tank 2']; 
x_label = 'Time [min]'; y_label = 'Height [cm]';
leg = ["10 %" ; "25 %" ; "50 %"];% Legends
col = ["r";"b";"k"];

for i = 1:2
    [T_10,X_10] = ode15s(@FourTankSystem,[t_0 t_end],x_s,[],F_in(1,:,i+1),p);
    [T_25,X_25] = ode15s(@FourTankSystem,[t_0 t_end],x_s,[],F_in(2,:,i+1),p);
    [T_50,X_50] = ode15s(@FourTankSystem,[t_0 t_end],x_s,[],F_in(3,:,i+1),p);
    
    % Save data to vectors
    T_10 = [T ; T_10]; T_25 = [T ; T_25]; T_50 = [T ; T_50];
    X_10 = [X ; X_10]; X_25 = [X ; X_25]; X_50 = [X ; X_50];

    % Generate plot
    % Save data in cells
    STATES = cell(1,3);
    STATES(1,1) = num2cell(X_10/(rho*A1),[1 2]);
    STATES(1,2) = num2cell(X_25/(rho*A1),[1 2]);
    STATES(1,3) = num2cell(X_50/(rho*A1),[1 2]);
    TIME = cell(1,3);
    TIME(1,1) = num2cell(T_10/60,[1 2]);
    TIME(1,2) = num2cell(T_25/60,[1 2]);
    TIME(1,3) = num2cell(T_50/60,[1 2]);
    sgTit = [];

    plotFunc(STATES,TIME,n_subplot,StateNames,x_label,y_label,leg,sgTit,col)
end % i 

%% Problem 3.2 Step Stochastic
% Measurement noise
L_y_low = chol([0.2 0 ; 0 0.4])';
L_y_med = chol([1 0 ; 0 2])';
L_y_hig = chol([5 0 ; 0 10])';

% State noise
% L_x_low = chol([diag([0.2 0.4 0.6 0.8])]);
% L_x_med = chol([diag([2 4 6 8])]);
% L_x_hig = chol([diag([20 40 60 80])]);

L_x_low = diag([0 0 5 10]);
L_x_med = diag([0 0 50 100]);
L_x_hig = diag([0 0 500 1000]);

t = t0:Ts:t_f;
x = zeros(4,length(t)); x(:,1) = x_s';
x = x_s';
y = zeros(2,length(t));

% First generate "steady state" stochastic simulation
for k = 1:t_f/Ts
    % Generate random noise on states and outputs
    x_LN = mvnrnd(zeros(size(A_c,1),1),L_x_low)';
    x_MN = mvnrnd(zeros(size(A_c,1),1),L_x_med)';
    x_HN = mvnrnd(zeros(size(A_c,1),1),L_x_hig)';
    y_LN = mvnrnd(zeros(size(C_c,1),1),L_y_low)';
    y_MN = mvnrnd(zeros(size(C_c,1),1),L_y_med)';
    y_HN = mvnrnd(zeros(size(C_c,1),1),L_y_hig)';
    
    %Determine the outputs with different noise levels
    y(:,k) = x(1:2,k)/(rho*A1);
    y0_L(:,k) = y(:,k) + y_LN;
    y0_M(:,k) = y(:,k) + y_MN;
    y0_H(:,k) = y(:,k) + y_HN;
    
    % Run simulation - Only use 1 value for the states
    [~,x_i] = ode15s(@FourTankSystem,[t(k) t(k+1)],x(:,k),[],F_in(1,:,1),p);
    x(:,k+1) = x_i(end,:)';

    % Determine states with different noise levels
    x0_L(:,k+1) = x(:,k+1) + x_LN;
    x0_M(:,k+1) = x(:,k+1) + x_MN;
    x0_H(:,k+1) = x(:,k+1) + x_HN;
end % k


n_subplot = 2; % Number of subplots
StateNames = ['Tank 1';'Tank 2']; 
x_label = 'Time [min]'; y_label = 'Height [cm]';
% n_subplot = 4; % Number of subplots
% StateNames = ['Tank 3'; 'Tank 4';'Tank 1';'Tank 2']; 
% x_label = 'Time [min]'; y_label = 'mass [g]';
leg = ["High" ; "Medium" ; "Low"];% Legends
col = ["r";"b";"k"];
for i = 1:3
    for k = 1:t_f/Ts
        % Generate random noise on states and outputs
        x_LN = mvnrnd(zeros(size(A_c,1),1),L_x_low)';
        x_MN = mvnrnd(zeros(size(A_c,1),1),L_x_med)';
        x_HN = mvnrnd(zeros(size(A_c,1),1),L_x_hig)';
        y_LN = mvnrnd(zeros(size(C_c,1),1),L_y_low)';
        y_MN = mvnrnd(zeros(size(C_c,1),1),L_y_med)';
        y_HN = mvnrnd(zeros(size(C_c,1),1),L_y_hig)';
        
        %Determine the outputs with different noise levels
        y(:,k) = x(1:2,k)/(rho*A1);
        y_L(:,k) = y(:,k) + y_LN;
        y_M(:,k) = y(:,k) + y_MN;
        y_H(:,k) = y(:,k) + y_HN;
        
        % Run simulation - Only use 1 value for the states
        [~,x_i] = ode15s(@FourTankSystem,[t(k) t(k+1)],x(:,k),[],F_in(i,:,2),p);
        x(:,k+1) = x_i(end,:)';

        % Determine states with different noise levels
        x_L(:,k+1) = x(:,k+1) + x_LN;
        x_M(:,k+1) = x(:,k+1) + x_MN;
        x_H(:,k+1) = x(:,k+1) + x_HN;
    end % k
    
    % used in norm step
    if i == 1
        y_norm_step(:,[1,4]) = y_L';
    end
    if i == 2
        y_norm_step(:,[2,5]) = y_L';
    end
    if i == 3
        y_norm_step(:,[3,6]) = y_L';
    end
   
    
    % Save data to vectors
    t1 = (t0:Ts:t_f-1)/60; % In order to have the correct size
    t2 = (t_f:Ts:2*t_f-1)/60; % In order to have the correct size
    T = [t1' ; t2'];
    Y_L = [y0_L' ; y_L'];
    Y_M = [y0_M' ; y_M'];
    Y_H = [y0_H' ; y_H'];

%     X_L = [x0_L(:,1:end-1)' ; x_L(:,1:end-1)'];
%     X_M = [x0_M(:,1:end-1)' ; x_M(:,1:end-1)'];
%     X_H = [x0_H(:,1:end-1)' ; x_H(:,1:end-1)'];

    % Generate plot
    % Save data in cells
    STATES = {Y_H , Y_M , Y_L};
%     STATES = {X_H , X_M , X_L};
    TIME = {T , T , T};
    sgTit = [];

    plotFunc(STATES,TIME,n_subplot,StateNames,x_label,y_label,leg,sgTit,col)
end % i 
   

%% Problem 3.3 Normalized step
% Y1_stoc =  [Y_L(end/2+1:end,1)  Y_M(end/2+1:end,1) Y_H(end/2+1:end,1)];
% Y2_stoc =  [Y_L(end/2+1:end,2)  Y_M(end/2+1:end,2) Y_H(end/2+1:end,2)];

for i = 1:3 % Only look at input 1
    [T,X] = ode15s(@FourTankSystem,[t0:1:t_f],x_s,[],F_in(i,:,2),p);
    Step_y1 = (X(:,1)-x_s(1))/(F_in(i,1,2)-F_in(1,1,1));
    Step_y2 = (X(:,2)-x_s(2))/(F_in(i,1,2)-F_in(1,1,1)); 
    Y1_stoc_norm(:,i) = ( y_norm_step(:,i) - (x_s(1) / (rho*A1)) )/(F_in(i,1,2)-F_in(1,1,1)); 
    Y2_stoc_norm(:,i) = ( y_norm_step(:,3+i) - (x_s(2) / (rho*A1)) )/(F_in(i,1,2)-F_in(1,1,1)); 
    step_u1(i,:) = {T , Step_y1 , Step_y2};
end

for i = 1:3 % Only look at input 2
    [T,X] = ode15s(@FourTankSystem,[t0:1:t_f],x_s,[],F_in(i,:,3),p);
    Step_y1 = (X(:,1)-x_s(1))/(F_in(i,2,3)-F_in(i,2,1));
    Step_y2 = (X(:,2)-x_s(2))/(F_in(i,2,3)-F_in(i,2,1));
    step_u2(i,:) = {T , Step_y1 , Step_y2};
end

figure
subplot(2,2,1); hold on
plot(step_u1{1,1}/60,step_u1{1,2}/(rho*A1),'r')
plot(step_u1{2,1}/60,step_u1{2,2}/(rho*A1),'b')
plot(step_u1{3,1}/60,step_u1{3,2}/(rho*A1),'k')
legend('10 %','25 %','50 %')
hold off; grid on; ylabel('$h_1$'); title('$u_1$'); 

subplot(2,2,2); hold on
plot(step_u2{1,1}/60,step_u2{1,2}/(rho*A1),'r')
plot(step_u2{2,1}/60,step_u2{2,2}/(rho*A1),'b')
plot(step_u2{3,1}/60,step_u2{3,2}/(rho*A1),'k')
hold off; grid on; title('$u_2$');

subplot(2,2,3); hold on
plot(step_u1{1,1}/60,step_u1{1,3}/(rho*A1),'r')
plot(step_u1{2,1}/60,step_u1{2,3}/(rho*A1),'b')
plot(step_u1{3,1}/60,step_u1{3,3}/(rho*A1),'k')
hold off; grid on; ylabel('$h_2$'); xlabel('Time [min]');

subplot(2,2,4); hold on
plot(step_u2{1,1}/60,step_u2{1,3}/(rho*A1),'r')
plot(step_u2{2,1}/60,step_u2{2,3}/(rho*A1),'b')
plot(step_u2{3,1}/60,step_u2{3,3}/(rho*A1),'k')
hold off; grid on; xlabel('Time [min]');

% Stochastics normalized!
figure
subplot(1,2,1); hold on
plot(t1',Y1_stoc_norm(:,1),'r')
plot(t1',Y1_stoc_norm(:,2),'b')
plot(t1',Y1_stoc_norm(:,3),'k')
legend('10 %','25 %','50 %')
hold off; grid on; ylabel('$h_1$'); title('$u_1$'); xlabel('Time [min]');

subplot(1,2,2); hold on
plot(t1',Y2_stoc_norm(:,1),'r')
plot(t1',Y2_stoc_norm(:,2),'b')
plot(t1',Y2_stoc_norm(:,3),'k')
hold off; grid on; ylabel('$h_2$'); xlabel('Time [min]'); title('$u_1$'); 

%% Problem 3.4 - Identify transfer functions
% u1 y1
num_u1y1 = step_u1{1,2}(end)/(rho*A1);
den_u1y1 = [2.3*60 1]; % for num*0.63
theta = 0;
[s_u1y1,ts_u1y1]=sisoctf2dstep(num_u1y1,den_u1y1,theta,1,t_f);

% u2 y1
num_u2y1 = step_u2{1,2}(end)/(rho*A1);
p1 = 90;
p2 = 138;
den_u2y1 = [p1*p2 p1+p2 1]; 
theta = 0;
[s_u2y1,ts_u2y1]=sisoctf2dstep(num_u2y1,den_u2y1,theta,1,t_f);

% u1 y2
num_u1y2 = step_u1{1,3}(end)/(rho*A1);
p1 = 80;
p2 = 173;
den_u1y2 = [p1*p2 p1+p2 1]; 
theta = 0;
[s_u1y2,ts_u1y2]=sisoctf2dstep(num_u1y2,den_u1y2,theta,1,t_f);

% u2 y2
num_u2y2 = step_u2{1,3}(end)/(rho*A1);
den_u2y2 = [2.55*60 1]; % for num*0.63
theta = 0;
[s_u2y2,ts_u2y2]=sisoctf2dstep(num_u2y2,den_u2y2,theta,1,t_f);

figure
subplot(2,2,1); hold on
plot(ts_u1y1/60,s_u1y1,'b'); plot(step_u1{1,1}/60,step_u1{1,2}/(rho*A1),'--r')
hold off; grid on; ylabel('$h_1$'); title('$u_1$'); legend('Fitted','Real','Location','best')
subplot(2,2,2); hold on
plot(ts_u2y1/60,s_u2y1,'b'); plot(step_u2{1,1}/60,step_u2{1,2}/(rho*A1),'--r')
hold off; grid on; title('$u_2$');
subplot(2,2,3); hold on
plot(ts_u1y2/60,s_u1y2,'b'); plot(step_u1{1,1}/60,step_u1{1,3}/(rho*A1),'--r')
hold off; grid on; ylabel('$h_2$');
subplot(2,2,4); hold on; xlabel('Time [sec]')
plot(ts_u2y2/60,s_u2y2,'b'); plot(step_u2{1,1}/60,step_u2{1,3}/(rho*A1),'--r')
hold off; grid on; xlabel('Time [sec]')

%% Problem 3.5 Accuracy
% Determine precision of estimate TF
e_u1y1 = step_u1{1,2}/(rho*A1)-s_u1y1; avg_u1y1 = mean(e_u1y1), var_u1y1 = var(e_u1y1), per_u1y1 = avg_u1y1/s_u1y1(end)*100
e_u2y1 = step_u2{1,2}/(rho*A1)-s_u2y1; avg_u2y1 = mean(e_u2y1), var_u2y1 = var(e_u2y1), per_u2y1 = avg_u2y1/s_u2y1(end)*100
e_u1y2 = step_u1{1,3}/(rho*A1)-s_u1y2; avg_u1y2 = mean(e_u1y2), var_u1y2 = var(e_u1y2), per_u1y2 = avg_u1y2/s_u1y2(end)*100
e_u2y2 = step_u2{1,3}/(rho*A1)-s_u2y2; avg_u2y2 = mean(e_u2y2), var_u2y2 = var(e_u2y2), per_u2y2 = avg_u2y2/s_u2y2(end)*100

%% Problem 3.6 Impulse Response Coefficients
% Look at sec 1.9 in "LinearModelPredictiveControlToolbox"
% Observer canonical form based on normalized TF
tf_11 = tf(num_u1y1,den_u1y1);
temp = den_u1y1(1); % To normalize
Aoc_11 = -den_u1y1(2) / temp;
Boc_11 = num_u1y1 / temp;
Coc_11 = 1;
Doc_11 = 0;

tf_12 = tf(num_u1y2,den_u1y2);
temp = den_u1y2(1); % To normalize
Aoc_12 = [-den_u1y2(2)/temp     1;
                 -den_u1y2(3)/temp      0];
Boc_12 = [0 ; 
                  num_u1y2/temp];
Coc_12 =  [1                                0];
Doc_12 = 0;

tf_21 = tf(num_u2y1,den_u2y1);
temp = den_u2y1(1); % To normalize
Aoc_21 = [-den_u2y1(2)/temp     1;
                 -den_u2y1(3)/temp      0];
Boc_21 = [0 ; 
                  num_u2y1/temp];
Coc_21 =  [1                                0];
Doc_21 = 0;

tf_22 = tf(num_u2y2,den_u2y2);
temp = den_u2y2(1); % To normalize
Aoc_22 = -den_u2y2(2) / temp;
Boc_22 = num_u2y2 / temp;
Coc_22 = 1;
Doc_22 = 0;

% Discretize
Ts = 15; % Given in assignment
Q = expm([Aoc_11 Boc_11 ; 0 0]*Ts ); A_11 = Q(1,1) ; B_11 = Q(1,2);
Q = expm([Aoc_12 Boc_12 ; 0 0 0]*Ts ); A_12 = Q(1:2,1:2) ; B_12 = Q(1:2,3);
Q = expm([Aoc_21 Boc_21 ; 0 0 0]*Ts ); A_21 = Q(1:2,1:2) ; B_21 = Q(1:2,3);
Q = expm([Aoc_22 Boc_22 ; 0 0]*Ts ); A_22 = Q(1,1) ; B_22 = Q(1,2);

% Determine Hankel matrix
H0 = [Doc_11 , Doc_12 ; Doc_21 , Doc_22];
for i = 0:6
    pil = i*2;
    H(:,pil+1:pil+2) = [Coc_11*A_11^i*B_11 , Coc_12*A_12^i*B_12 ; Coc_21*A_21^i*B_21 , Coc_22*A_22^i*B_22];
end % i

% Determine which hankel norm ensuring minimal representation
Hankel2 = [H(:,1:4); H(:,3:6)]; rankH2 = rank(Hankel2)
Hankel3 = [H(:,1:6) ; H(:,3:8) ; H(:,5:10)]; rankH3 = rank(Hankel3)
Hankel4 = [H(:,1:8) ; H(:,3:10) ; H(:,5:12) ; H(:,7:14)]; rankH4 = rank(Hankel4)

% The rank is the same! Hankel3 is a minimal representation!
[K,lambda,L] = svd(Hankel3);
Kappa_1 = K(:,1:4)
lambda_1 = lambda(1:4,1:4)
L_1 = L(:,1:4)
Kappa_1*lambda_1*L_1'

%New MIMO system (look from H_2 to H_2*N-1)
H_tilde = [H(:,5:10) ; H(:,7:12) ; H(:,9:14)]

% Now determine state space contributions
A_k_tf = inv(sqrtm(lambda_1))*Kappa_1'*H_tilde*L_1*inv(sqrtm(lambda_1))
B_k_tf = sqrtm(lambda_1)*(L_1(1:2,:))'
C_k_tf = Kappa_1(1:2,:)*sqrtm(lambda_1)
D_k_tf = H0


N = t_f/Ts; % Number of samples
t = 0:Ts:t_f;
[x11 , x12 , x21 , x22] = MarkovParamters(A_k_tf,B_k_tf,C_k_tf,D_k_tf,N);
figure
hold on
plot(t/60,x11,'r'); plot(t/60,x12,'b'); plot(t/60,x21,'m'); plot(t/60,x22,'g')
hold off; grid on
legend('u_1y_1','u_1y_2','u_2y_1','u_2y_2')
xlabel('Time [sec]'); ylabel('Height [cm]')

%% Problem 4.1 Linearize
% Look at "FourTankSys_SS"
A_alg,B_alg,C_alg,E_alg, D_alg % Algebraic matrix
x_s % steady state
A_c, B_c, C_c,E_c,D_c % Numeric matrix

%% Problem 4.2 Poles, Gains and Zeros
% Poles
lambda = eig(A_c)

% gains
K_11 = -C_c(1,:)*inv(A_c)*B_c(:,1)
K_12 = -C_c(1,:)*inv(A_c)*B_c(:,2)
K_21 = -C_c(2,:)*inv(A_c)*B_c(:,1)
K_22 = -C_c(2,:)*inv(A_c)*B_c(:,2)


% Zeros
N = [eye(size(A_c)) , zeros(size(A_c,1),1) ; zeros(1,5)];

M_11 = [A_c , B_c(:,1) ; C_c(1,:) , zeros(1,1)]; z_11 = eig(M_11,N);% All zeros which is infinite should be removed 
z_11 = z_11(isfinite(z_11)) % Continious zeros

M_12 = [A_c , B_c(:,2) ; C_c(1,:) , zeros(1,1)]; z_12 = eig(M_12,N);% All zeros which is infinite should be removed 
z_12 = z_12(isfinite(z_12)) % Continious zeros

M_21 = [A_c , B_c(:,1) ; C_c(2,:) , zeros(1,1)]; z_21 = eig(M_21,N);% All zeros which is infinite should be removed 
z_21 = z_21(isfinite(z_21)) % Continious zeros

M_22 = [A_c , B_c(:,2) ; C_c(2,:) , zeros(1,1)]; z_22 = eig(M_22,N);% All zeros which is infinite should be removed 
z_22 = z_22(isfinite(z_22)) % Continious zeros

%% Problem 4.3 Transfer functions
s = tf('s');
TF_lin = C_c*(s*eye(4)-A_c)^-1*B_c

%% Problem 4.4 Comparsion
% Time constant linear model
tau_u1y1 = abs(1./eig(TF_lin(1,1)))
tau_u2y1 = abs(1./eig(TF_lin(1,2)))
tau_u2y1 = abs(1./eig(TF_lin(2,1)))
tau_u2y2 = abs(1./eig(TF_lin(2,2)))

%% Problem 4.5 Discretaization
% Look at "FourTankSys_SS"
A_k,B_k,E_k

%% Problem 4.6 Markov - linear and step comparison
N = t_f/Ts; % Number of samples
t = 0:Ts:t_f;
[x11_tf , x12_tf , x21_tf , x22_tf] = MarkovParamters(A_k_tf,B_k_tf,C_k_tf,D_k_tf,N); % Response from TF determine in pr 3.6
[x11_ss , x12_ss , x21_ss , x22_ss] = MarkovParamters(A_k_tf,B_k_tf,C_k_tf,D_k_tf,N); % Response from SS determine in pr 4.5


figure
hold on
plot(t/60,x11_ss,'r'); plot(t/60,x11_tf,'--k'); 
plot(t/60,x12_ss,'b'); plot(t/60,x12_tf,'--k'); 
plot(t/60,x21_ss,'m'); plot(t/60,x21_tf,'--k'); 
plot(t/60,x22_ss,'g'); plot(t/60,x22_tf,'--k'); 
hold off; grid on
legend('u_1y_1','','u_1y_2','','u_2y_1','','u_2y_2','TF')
xlabel('Time [sec]'); ylabel('Height [cm]')


%% Problem 5.2 Static and Dynamic Kalman Filter 
h_s = (x_s(:,1:2)/(rho*A1))';
sys = ss(A_k,[B_k,E_k],C_c,[D_c,zeros(size(E_c,2))]);
mat = expm([-A_c E_c*E_c' ; zeros(size(A_c)) A_c']*Ts);
phi_12 = mat(1:4,5:8);
phi_22 = mat(5:8,5:8);
Q_hat = phi_22'*phi_12;

[sys_aug,Q_hat_aug] = sys_aug_func(sys,Q_hat);

% noise
sigma_R = 0.5; % COVARIANCE
Rd = 5;
v_k = mvnrnd(zeros(t_f/Ts,size(C_c,1)),eye(2)*sigma_R)';
d_k = mvnrnd(zeros(t_f/Ts,size(E_c,2)),eye(2)*Rd)';

u_kal = u*0.1;
% STATIC
[x,y,x_hat_sta,y_hat_sta] = KalmanSta(sys,sys_aug,Q_hat,Q_hat_aug,sigma_R,Rd,t,x0',v_k,d_k,u_kal);

% DYNAMIC
[~,~,x_hat_dyn,y_hat_dyn] = KalmanDyn(sys,sys_aug,Q_hat,Q_hat_aug,sigma_R,Rd,t,x0',v_k,d_k,u_kal);

t_plot = (1:Ts:t_f)/60;

STATES = {(x+x_s')', (x_hat_sta+[x_s' ; 0 ; 0])', (x_hat_dyn+[x_s' ; 0 ; 0])'}; % Add more arrays and split by COMMA
TIME = {t_plot , t_plot , t_plot}; % Add more arrays and split by COMMA
sgTit = []; 
leg = ["Measured" ; "Static Kalman" ; "Dynamic Kalman"]; % Legends
x_label = 'Time [min]'; y_label = 'Mass [g]';
n_subplot = 4;
StateNames = ['Tank 3'; 'Tank 4';'Tank 1';'Tank 2']; 
col = ["k";"r";"--b"];
plotFunc(STATES,TIME,n_subplot,StateNames,x_label,y_label,leg,sgTit,col)

STATES = {(y+h_s)' , (y_hat_sta+h_s)' , (y_hat_dyn+h_s)'}; % Add more arrays and split by COMMA
n_subplot = 2;
StateNames = ['Tank 1';'Tank 2'];
y_label = 'Height [cm]';
plotFunc(STATES,TIME,n_subplot,StateNames,x_label,y_label,leg,sgTit,col)


%% Problem 5.3 Static and Dynamic Kalman Filter - STEP CHANGES
% Run initla steady state simulation
% STATIC
d_k = mvnrnd(zeros(t_f/Ts,size(E_c,2)),eye(2)*Rd)';
[x,y,x_hat_sta,y_hat_sta] = KalmanSta(sys,sys_aug,Q_hat,Q_hat_aug,sigma_R,Rd,t,x0',v_k,d_k,u_kal);
% DYNAMIC
[~,~,x_hat_dyn,y_hat_dyn] = KalmanDyn(sys,sys_aug,Q_hat,Q_hat_aug,sigma_R,Rd,t,x0',v_k,d_k,u_kal);
x_s_kal = x(:,end)'; % Set final value as new steady state

d_k = 250*0.1 + mvnrnd(zeros(t_f/Ts,size(E_c,2)),eye(2)*Rd)';
% STATIC
[x_10,y_10,x_hat_sta_10,y_hat_sta_10] = KalmanSta(sys,sys_aug,Q_hat,Q_hat_aug,sigma_R,Rd,t,x_s_kal,v_k,d_k,u_kal);
% DYNAMIC
[~,~,x_hat_dyn_10,y_hat_dyn_10] = KalmanDyn(sys,sys_aug,Q_hat,Q_hat_aug,sigma_R,Rd,t,x_s_kal,v_k,d_k,u_kal);

% Save data to vectors
% T_10 = [t0:Ts:t_f , t_f:Ts:t_f*2]; 
x_10 = [x , x_10]; % Measured
x_hat_sta_10 = [x_hat_sta , x_hat_sta_10];  % Static kalman
x_hat_dyn_10 = [x_hat_dyn , x_hat_dyn_10]; % dynamic kalman

y_10 = [y , y_10];
y_hat_sta_10 = [y_hat_sta , y_hat_sta_10];
y_hat_dyn_10 = [y_hat_dyn , y_hat_dyn_10];

t_plot = (1:Ts:t_f*2)/60;

STATES = {(x_10+x_s')', (x_hat_sta_10+[x_s' ; 0 ; 0])', (x_hat_dyn_10+[x_s' ; 0 ; 0])'}; % Add more arrays and split by COMMA
TIME = {t_plot , t_plot , t_plot}; % Add more arrays and split by COMMA
n_subplot = 4; % Number of subplots
StateNames = ['Tank 3'; 'Tank 4';'Tank 1';'Tank 2']; 
x_label = 'Time [min]'; y_label = 'Mass [g]';
leg = ["Measured" ; "Static Kalman" ; "Dynamic Kalman"]; % Legends
col = ["k";"r";"--b"];
sgTit = [];
plotFunc(STATES,TIME,n_subplot,StateNames,x_label,y_label,leg,sgTit,col)

STATES = {(y_10+h_s)' , (y_hat_sta_10+h_s)' , (y_hat_dyn_10+h_s)'}; % Add more arrays and split by COMMA
n_subplot = 2;
StateNames = ['Tank 1';'Tank 2'];
y_label = 'Height [cm]';
plotFunc(STATES,TIME,n_subplot,StateNames,x_label,y_label,leg,sgTit,col)

%% Problem 5.4 Kalman non-linear
A_aug = sys_aug.A; C_aug = sys_aug.C;

% Static
R = eye(size(C_c,1),size(C_c,1)) * sqrt(sigma_R);
P = idare(A_aug',C_aug',Q_hat_aug,R);
P_sta = P;
Re = C_aug*P*C_aug'+R;
L = P*C_aug'*inv(Re); % Filter can also use DLQE

% Data vector!    
x_now = x_s'; x_hat_sta = [x0 ; zeros(2,1)]; x_hat_dyn = x_hat_sta;
y = h_s; y_hat_sta = y; y_hat_dyn = y;
data=[0 x_now' y' x_hat_sta'+[x_s 0 0] y_hat_sta' x_hat_dyn'+[x_s 0 0] y_hat_dyn'];

for i = 1:Ts:t_f
    v_k = mvnrnd(zeros(size(C_c,1),1),sigma_R);
    d_k = mvnrnd(zeros(size(E_c,2),1),Rd);
    F_dev = [u_kal' , d_k']; F_stoc = [u' d'] + F_dev;
    
    % Non linear sim
    t = []; x_nl = [];
    [t,x_nl] = ode15s(@FourTankSystem,[i i+Ts],x_now,[],F_stoc,p);

    % Update including noise addition
    x_now = x_nl(end,:)';
    y = x_now(1:2,:)/(rho*A1) + v_k;

    % Kalman
    y_dev = y - h_s; % Kalman needs deviation variables
    type = 1; % 1 for static 2 for dynamic
    [x_hat_sta,y_hat_sta,P_sta] = oneStepKalman(sys_aug,Q_hat_aug,sigma_R,P_sta,R,Re,L,x_hat_sta,y_dev,F_dev(1:2),type);
    type = 2; % 1 for static 2 for dynamic
    [x_hat_dyn,y_hat_dyn,P] = oneStepKalman(sys_aug,Q_hat_aug,sigma_R,P,R,Re,L,x_hat_sta,y_dev,F_dev(1:2),type);
    data = [data; t(end,:) x_now' y' x_hat_sta'+[x_s 0 0] y_hat_sta'+h_s' x_hat_dyn'+[x_s 0 0] y_hat_dyn'+h_s'];
end %i

% NOW STEP ON DISTURBANCE!
for i = 1:Ts:t_f
    v_k = mvnrnd(zeros(size(C_c,1),1),sigma_R);
    d_k = mvnrnd(zeros(size(E_c,2),1),Rd);
    F_dev = [u_kal' , d_k']; F_stoc = [u' d'*1.1] + F_dev;
    
    % Non linear sim
    t = []; x_nl = [];
    [t,x_nl] = ode15s(@FourTankSystem,[i i+Ts],x_now,[],F_stoc,p);

    % Update including noise addition
    x_now = x_nl(end,:)';
    y = x_now(1:2,:)/(rho*A1) + v_k;

    % Kalman
    y_dev = y - h_s; % Kalman needs deviation variables
    type = 1; % 1 for static 2 for dynamic
    [x_hat_sta,y_hat_sta,P_sta] = oneStepKalman(sys_aug,Q_hat_aug,sigma_R,P_sta,R,Re,L,x_hat_sta,y_dev,F_dev(1:2),type);
    type = 2; % 1 for static 2 for dynamic
    [x_hat_dyn,y_hat_dyn,P] = oneStepKalman(sys_aug,Q_hat_aug,sigma_R,P,R,Re,L,x_hat_sta,y_dev,F_dev(1:2),type);
    data = [data; t(end,:) x_now' y' x_hat_sta'+[x_s 0 0] y_hat_sta'+h_s' x_hat_dyn'+[x_s 0 0] y_hat_dyn'+h_s'];
end %i
t_plot = (1:Ts:t_f*2)/60;

STATES ={data(1:end-1,2:5) , data(1:end-1,8:11) , data(1:end-1,16:19)};
TIME = {t_plot , t_plot , t_plot}; % Add more arrays and split by COMMA
n_subplot = 4; % Number of subplots
StateNames = ['Tank 3'; 'Tank 4';'Tank 1';'Tank 2']; 
x_label = 'Time [min]'; y_label = 'Mass [g]';
leg = ["Measured" ; "Static Kalman" ; "Dynamic Kalman"]; % Legends
col = ["k";"r";"--b"];
sgTit = [];
plotFunc(STATES,TIME,n_subplot,StateNames,x_label,y_label,leg,sgTit,col)

STATES ={data(1:end-1,6:7) , data(1:end-1,14:15) , data(1:end-1,22:23)};
n_subplot = 2;
StateNames = ['Tank 1';'Tank 2'];
y_label = 'Height [cm]';
plotFunc(STATES,TIME,n_subplot,StateNames,x_label,y_label,leg,sgTit,col)

%% Problem 7 Unconstrained MPC
Q = 100;% Tuning parameter
S = eye(2)*0.01;% Tuning parameter
t_R = 10*60; % The reference is constant for 10 minutes
t_f = t_R*3; % New time in step reference
N = t_R/Ts; % Should be at the same length at time horizon of the dynamics

MPC_sys = MPCDesign(sys,Q,S,N);

v_k = mvnrnd(zeros(t_f/Ts,size(C_c,1)),eye(2)*sigma_R)';
d_k = mvnrnd(zeros(t_f/Ts,size(E_c,2)),eye(2)*Rd)';

x = zeros(4,t_f/Ts+1);
u_vec = zeros(2,t_f/Ts+1);
R1 = [0 ; 0];
R2 = [10 ; 10];
R3 = [20 ;-10];

R=zeros(2*(N+t_f/Ts),1);
R(1:2*(t_f/Ts/3),1)=kron(ones(t_f/Ts/3,1),R1);
R(2*(t_f/Ts/3)+1:2*(2*t_f/Ts/3),1)=kron(ones(t_f/Ts/3,1),R2);
R(2*(2*t_f/Ts/3)+1:2*(N+t_f/Ts),1)=kron(ones(N+t_f/Ts/3,1),R3);

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

%% Problem 8 Input constrained MPC
ub_c = F_in(1,1,1)-150; ub = ub_c*ones(2*N,1);
lb_c = -F_in(1,1,1); lb = lb_c*ones(2*N,1);
u_rate = 10*ones(2*N,1);   % maximum rate  of change
l_rate = -10*ones(2*N,1);  % minimum rate of change

Q = 100;% Tuning parameter
S = eye(2)*0.001; % Tuning parameter

MPC_sys = MPCDesign(sys,Q,S,N);
x = zeros(4,t_f/Ts+1); 
u_vec = zeros(2,t_f/Ts+1);
cons = {lb,ub,l_rate,u_rate}; % Constraints

inputs = {x,x_s,u_vec,R,v_k,d_k};

type = 2; % Input constrained MPC
[y,y_nl,u_vec,u_nl] = MPC_Sim(sys,sys_aug,MPC_sys,[],Q_hat,Q_hat_aug,sigma_R,t_f,t_R,Ts,inputs,N,cons,[],[],p,type);

figure
subplot (2,1,1); hold on
plot(t_plot,y(1,:)+h_s(1),'r');plot(t_plot,y_nl(1,:),'--m'); plot(t_plot,R_plot(:,1),'--k'); plot(t_plot,y(2,:)+h_s(2),'b');plot(t_plot,y_nl(2,:),'--g'); plot(t_plot,R_plot(:,2),'-.k')
hold off; grid on; legend('h_1 - lin','h_1 - Nonlinear','r_1','h_2 - lin','h_2 - Nonlinear','r_2','Location','best'); ylabel('Height $[cm]$','Interpreter','latex')
subplot (2,1,2); hold on
plot(t_plot,u_vec(1,:)+u(1),'r');plot(t_plot,u_nl(1,:)+u(1),'--m'); plot(t_plot,u_vec(2,:)+u(2),'b');plot(t_plot,u_nl(2,:)+u(2),'--g'); yline(ub_c+300,'--k','LineWidth',2); yline(lb_c+300,'--k','LineWidth',2);
hold off; grid on; legend('u_1','u_1 - Nonlinear','u_2','u_2 - Nonlinear','Location','best'); xlabel('Time [min]'), ylabel('Flow $[cm^3/s]$'); ylim([-100 700])


%% Problem 9 Hard input soft output constrained MPC
W_z = 100*eye(2);
W_u = 0*eye(2);
W_du = 0*eye(2);
W = [W_z W_u W_du];

% Boundaries
W_t1 = 10000*eye(2);
W_t2 = 10000*eye(2);
W_s1 = 10000*eye(2);
W_s2 = 10000*eye(2);
soft_cons = {W_t1,W_t2,W_s1,W_s2};

MPC_sys_soft = Soft_MPCDesign(sys,MPC_sys,N,W,soft_cons);

% Boundaries Inputs
U_bar = kron(ones(N,1),[0;0]);
ub_c = F_in(1,1,1)-150; U_max = ub_c*ones(2*N,1);
lb_c = -F_in(1,1,1); U_min = lb_c*ones(2*N,1);

% Boundaries input deviation
D_u_min = kron(ones(N,1),[-10;-10]);
D_u_max = kron(ones(N,1),[10;10]);

% Boundaries output
Z_max = kron(ones(N,1),[110-h_s(1);120-h_s(2)]);
Z_min = kron(ones(N,1),[85-h_s(1);105-h_s(2)]);

cons_InOut = {U_bar,U_min,U_max,D_u_min,D_u_max,Z_max,Z_min};

x = zeros(4,t_f/Ts+1); 
u_vec = zeros(2,t_f/Ts+1);
inputs = {x,x_s,u_vec,R,v_k,d_k};
type = 3; % Input/Output constrained MPC
[y,y_nl,u_vec,u_nl] = MPC_Sim(sys,sys_aug,MPC_sys,MPC_sys_soft,Q_hat,Q_hat_aug,sigma_R,t_f,t_R,Ts,inputs,N,cons,cons_InOut,[],p,type);

figure
subplot (2,2,1); hold on
plot(t_plot,y(1,:)+h_s(1),'r');plot(t_plot,y_nl(1,:),'--b'); plot(t_plot,R_plot(:,1),'--k'); yline(Z_max(1,1)+h_s(1),'--m','LineWidth',2); yline(Z_min(1,1)+h_s(1),'--m','LineWidth',2);
hold off; grid on; legend('Linear','Nonlinear','r','Bounds','Location','best','Interpreter','latex'); ylabel('Height $[cm]$','Interpreter','latex'); ylim([80 130])
subplot(2,2,2); hold on
plot(t_plot,y(2,:)+h_s(2),'r');plot(t_plot,y_nl(2,:),'--b'); plot(t_plot,R_plot(:,2),'-.k'); yline(Z_max(2,1)+h_s(2),'--m','LineWidth',2); yline(Z_min(2,1)+h_s(2),'--m','LineWidth',2);
hold off; grid on; ylabel('Height $[cm]$','Interpreter','latex');ylim([80 130])
subplot (2,2,3); hold on
plot(t_plot,u_vec(1,:)+u(1),'r');plot(t_plot,u_nl(1,:)+u(1),'--b'); yline(ub_c+300,'--k','LineWidth',2); yline(lb_c+300,'--k','LineWidth',2);yline(ub_c+300,'--k','LineWidth',2); yline(lb_c+300,'--k','LineWidth',2);
hold off; grid on; legend('Linear','Nonlinear','Bounds','Location','best','Interpreter','latex'); xlabel('Time [min]'), ylabel('Flow $[cm^3/s]$'); ylim([-100 700])
subplot (2,2,4); hold on
plot(t_plot,u_vec(2,:)+u(2),'r');plot(t_plot,u_nl(2,:)+u(2),'--b'); yline(ub_c+300,'--k','LineWidth',2); yline(lb_c+300,'--k','LineWidth',2);yline(ub_c+300,'--k','LineWidth',2); yline(lb_c+300,'--k','LineWidth',2);
hold off; grid on;  xlabel('Time [min]'), ylabel('Flow $[cm^3/s]$'); ylim([-100 700])

%% Problem 11 Nonlinear MPC
t_f = 1800;
t = 0:Ts:t_f;
x = zeros(4,length(t)); x(:,1) = x_s';
x_hat = x;
x_hat1 = x;
u_cdekf = zeros(2,length(t)) + u; 
u_cdekf(1,40:80) = 400;  u_cdekf(1,80:121) = 500;
u_cdekf(2,40:80) = 200;  u_cdekf(2,80:121) = 250;
y = [];

v_k = mvnrnd(zeros(t_f/Ts+1,size(C_c,1)),eye(2)*sigma_R)'; % meas noise
d_k = mvnrnd(zeros(t_f/Ts+1,size(E_c,2)),eye(2)*Rd)'; % disturbance noise
P_k = eye(4);
Re=zeros(2,2,t_f/Ts+1);
for i = 1:t_f/Ts
    f = feval(@f_func,x(:,i),u_cdekf(:,i),d+d_k(:,i),p);
    g = feval(@sigma_func,x(:,i),p);
    x(:,i+1)=x(:,i)+f*Ts + g*v_k(:,i);
    y(:,i)=C_c*x(:,i) + v_k(:,i);

    [x_hat1(:,i),~,x_hat(:,i+1),P_k1,error(:,i),Re(:,:,i)] = CDEKF(E_k,@P_func,p,x_hat(:,i),u_cdekf(:,i),d+d_k(:,i),P_k,y(:,i),eye(2)*sigma_R,t,i,C_c);
    
end
y(:,t_f/Ts+1)=C_c*x(:,t_f/Ts+1) + v_k(:,t_f/Ts+1);

figure
subplot(2,2,[1 2]); hold on
% plot(t/60,x(1,:),'--k'); plot(t/60,x(2,:),'--k'); plot(t/60,x(3,:),'--k'); plot(t/60,x(4,:),'--k');
plot(t/60,x_hat(1,:),'r'); plot(t/60,x_hat(2,:),'b'); plot(t/60,x_hat(3,:),'g'); plot(t/60,x_hat(4,:),'m')
plot(t/60,x(1,:),'--k'); plot(t/60,x(2,:),'--k'); plot(t/60,x(3,:),'--k'); plot(t/60,x(4,:),'--k');
grid on; legend({'$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$','True states','','',''},'interpreter','latex'); ylabel('Mass $[g]$'); xlabel('Time $[min]$')
subplot(2,2,3)
plot(t/60,u_cdekf(:,:))
grid on; ylim([0 500]);legend('u_1','u_2'); ylabel('Flow $[cm^3/s]$'); xlabel('Time $[min]$')
subplot(2,2,4); hold on;
plot(t/60,x(1,:)-x_hat(1,:)); plot(t/60,x(2,:)-x_hat(2,:)); plot(t/60,x(3,:)-x_hat(3,:)); plot(t/60,x(4,:)-x_hat(4,:));
grid on; xlabel('Time $[min]$'); ylabel('Erorr $[g]$'); legend('e_1','e_2','e_3','e_4')

%% Problem 11.3 Parameter estimation
% import casadi.*
% % Estimation
%  %ipopt options
%  opts=struct;
%  opts.ipopt.max_iter=100;
%  opts.ipopt.print_level=0;
%  opts.print_time=0;
%  opts.ipopt.acceptable_tol=1e-8;
%  opts.ipopt.acceptable_obj_change_tol=1e-6;
% % Parameter optimization:
% A_1 = SX.sym('a1');
% A_2 = SX.sym('a2');
% w={};
% w0=[];
% lbw=[];
% ubw=[];
% % Continuous time dynamics
% Vm=0;
% for i=1:t_f/Ts
%     error=y(:,i)-[(1/A_1)*x_hat(1,i);(1/A_2)*x_hat(2,i)];
%     Vm=Vm+0.5*log(det(Re(:,:,i)))+0.5*error'*inv(Re(:,:,i))*error;
% end
% w = {w{:},A_1;A_2};
% lbw = [lbw; 0; 0;];
% ubw = [ubw;  inf;  inf;];
% w0 = [w0; 0;0;];
% 
% prob = struct('f', Vm, 'x', vertcat(w{:}));
% 
% solver = nlpsol('solver', 'ipopt', prob, opts);
% 
% % Solve the NLP
% sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw );
% estimated_param=full(sol.x)

%% Problem 12 Economic MPC
ub_c = F_in(1,1,1)-150; ub = ub_c*ones(2*N,1);
lb_c = -F_in(1,1,1); lb = lb_c*ones(2*N,1);
u_rate = 100*ones(2*N,1);   % maximum rate  of change
l_rate = -50*ones(2*N,1);  % minimum rate of change
Ubar=kron(ones(N,1),[0;0]);
cons = {lb,ub,l_rate,u_rate}; % Constraints

t_f = 1800; % 30 minutes
v_k = mvnrnd(zeros(t_f/Ts,size(C_c,1)),eye(2)*sigma_R)';
d_k = mvnrnd(zeros(t_f/Ts,size(E_c,2)),eye(2)*Rd)';
x = zeros(4,t_f/Ts+1);
u_vec = zeros(2,t_f/Ts+1);
R1 = [0 ; 0];
R2 = [10 ; 10];
R3 = [20 ;-10];
R=zeros(2*(N+t_f/Ts),1);
R(1:2*(t_f/Ts/3),1)=kron(ones(t_f/Ts/3,1),R1);
R(2*(t_f/Ts/3)+1:2*(2*t_f/Ts/3),1)=kron(ones(t_f/Ts/3,1),R2);
R(2*(2*t_f/Ts/3)+1:2*(N+t_f/Ts),1)=kron(ones(N+t_f/Ts/3,1),R3);
inputs = {x,x_s,u_vec,R,v_k,d_k};



c = [0.5 , 0.5]; rho_Econ = [20, 20]; % Large c makes u small. Large rho large makes V small
g_u=zeros(2*N,1); g_u(1:2:end) = c(1); g_u(2:2:end) = c(2);
g_v=zeros(2*N,1); g_v(1:2:end) = rho_Econ(1); g_v(2:2:end) = rho_Econ(2);



cons_Econ = [g_u g_v];

type = 4; % Economic
[y,y_nl,u_vec,u_nl] = MPC_Sim(sys,sys_aug,MPC_sys,[],Q_hat,Q_hat_aug,sigma_R,t_f,t_R,Ts,inputs,N,cons,[],cons_Econ,p,type);

t_plot = (1:Ts:t_f)/60;
R_plot = [ones(t_R/Ts,1)*R1(1)+h_s(1) , ones(t_R/Ts,1)*R1(2)+h_s(2) ;
               ones(t_R/Ts,1)*R2(1)+h_s(1) , ones(t_R/Ts,1)*R2(2)+h_s(2) ;
               ones(t_R/Ts,1)*R3(1)+h_s(1) , ones(t_R/Ts,1)*R3(2)+h_s(2)];


figure
subplot (2,1,1); hold on
plot(t_plot,y(1,:)+h_s(1),'r');plot(t_plot,y_nl(1,:),'--m'); plot(t_plot,R_plot(:,1),'--k'); plot(t_plot,y(2,:)+h_s(2),'b');plot(t_plot,y_nl(2,:),'--g'); plot(t_plot,R_plot(:,2),'-.k')
hold off; grid on; legend('$h_1$ - Linear','$h_1$ - Nonlinear','$r_1$','$h_2$ - Linear','$h_2$ - Nonlinear','$r_2$','Location','best','Interpreter','Latex'); ylabel('Height $[cm]$','Interpreter','latex')
subplot (2,1,2); hold on
plot(t_plot,u_vec(1,:)+u(1),'r');plot(t_plot,u_nl(1,:)+u(1),'--m'); plot(t_plot,u_vec(2,:)+u(2),'b');plot(t_plot,u_nl(2,:)+u(2),'--g'); 
hold off; grid on; legend('$u_1$ - Linear','$u_1$ - Nonlinear','$u_2$ - Linear','$u_2$ - Nonlinear','Location','best','Interpreter','Latex'); xlabel('Time [min]'), ylabel('Flow $[cm^3/s]$'); ylim([-100 700])



%% Problem  13 Classical control
t_f = 1800*2;
r_1 = zeros(t_f/Ts,1); r_1(end*1/3:end*2/3) = 10; r_1(end*2/3:end+1) = 20;
r_2 = zeros(t_f/Ts,1); r_2(end*1/3:end*2/3) = 10; r_2(end*2/3:end+1) =  -10;
r=[r_1';r_2'];

u_min= - 300;
u_max= 300;

x_P = zeros(4,t_f/Ts+1); x_PI = x_P; x_PID = x_P;
y_P = zeros(2,t_f/Ts+1); y_PI = y_P; y_PID = y_P;
u_P = zeros(2,t_f/Ts+1); u_PI = u_P; u_PID = u_P;
 
v_k = mvnrnd(zeros(t_f/Ts,size(C_c,1)),eye(2)*sigma_R)';
d_k = mvnrnd(zeros(t_f/Ts,size(E_c,2)),eye(2)*Rd)';

% P gains
K_P=[2 0;0 0.2];

% PI gains
K_PI = [7.5 0; 0 2.5]*0.1;
Ti = 2000;
i = [0;0];

% PID gains
 K_PID_P=[2 0;0 2];
 K_PID_I=[5 0; 0 5]*0.001;
 K_PID_D=[100 0; 0 100];
 Ti = 4000;
  i = [0;0];

for k=1:t_f/Ts
    if k == 1
        [u_P(:,k),e_P(:,k)]=P_con(r(:,k),y_P(:,k),u_P(:,k),u_min,u_max,K_P);
        [u_PI(:,k) e_PI(:,k)] = PI_con(r(:,k),y_PI(:,k),u_PI(:,k),u_min,u_max,K_PI,i,Ts,Ti);
    else
        [u_P(:,k),e_P(:,k)]=P_con(r(:,k),y_P(:,k-1),u_P(:,k-1),u_min,u_max,K_P);
        [u_PI(:,k) e_PI(:,k)] = PI_con(r(:,k),y_PI(:,k-1),u_PI(:,k-1),u_min,u_max,K_PI,i,Ts,Ti);
    end    
    
    if k <= 2
        [u_PID(:,k) e_PID(:,k)] = PID_con(r(:,k),u_PID(:,k),y_PID(:,k),y_PID(:,k),u_min,u_max,K_PID_P,i,K_PID_I,K_PID_D,Ts);
    else
        [u_PID(:,k) e_PID(:,k)] = PID_con(r(:,k),u_PID(:,k-1),y_PID(:,k-1),y_PID(:,k-2),u_min,u_max,K_PID_P,i,K_PID_I,K_PID_D,Ts);
    end

    x_P(:,k+1)=A_k*x_P(:,k)+B_k*u_P(:,k)+E_k*d_k(:,k);
    x_PI(:,k+1)=A_k*x_PI(:,k)+B_k*u_PI(:,k)+E_k*d_k(:,k);
    x_PID(:,k+1)=A_k*x_PID(:,k)+B_k*u_PID(:,k)+E_k*d_k(:,k);
    
    y_P(:,k)=C_c*x_P(:,k)+v_k(1:2,k);
    y_PI(:,k)=C_c*x_PI(:,k)+v_k(1:2,k);
    y_PID(:,k)=C_c*x_PID(:,k)+v_k(1:2,k);
    
end

[u_P(:,t_f/Ts+1),e_P(:,t_f/Ts+1)]=P_con(r(:,t_f/Ts+1),y_P(:,t_f/Ts),u_P(:,t_f/Ts),u_min,u_max,K_P);
[u_PI(:,t_f/Ts+1) e_PI(:,t_f/Ts+1)] = PI_con(r(:,t_f/Ts+1),y_PI(:,t_f/Ts),u_PI(:,t_f/Ts),u_min,u_max,K_PI,i,Ts,Ti);
[u_PID(:,t_f/Ts+1) e_PID(:,t_f/Ts+1)] = PID_con(r(:,t_f/Ts+1),u_PID(:,k-1),y_PID(:,t_f/Ts),y_PID(:,t_f/Ts),u_min,u_max,K_PID_P,i,K_PID_I,K_PID_D,Ts);

y_P(:,t_f/Ts+1)=C_c*x_P(:,t_f/Ts+1);
y_PI(:,t_f/Ts+1)=C_c*x_PI(:,t_f/Ts+1);
y_PID(:,t_f/Ts+1)=C_c*x_PID(:,t_f/Ts+1);

t_plot = (1:Ts:t_f+1)/60;

figure
subplot(2,3,1); hold on
plot(t_plot,y_P(1,:)+h_s(1),'r'); plot(t_plot,r(1,:)+h_s(1),'--k');
plot(t_plot,y_P(2,:)+h_s(2),'b'); plot(t_plot,r(2,:)+h_s(2),'-.k')
hold off; grid on; legend('h_1','r_1','h_2','r_2'); ylabel('Outputs'); title('P Controller'); ylim([80 140])
subplot(2,3,4); hold on
plot(t_plot,u_P(1,:)+300,'r'); plot(t_plot,u_P(2,:)+300,'b');
grid on; legend('u_1','u_2'); ylabel('Inputs'); ylim([0 600])

subplot(2,3,2); hold on
plot(t_plot,y_PI(1,:)+h_s(1),'r'); plot(t_plot,r(1,:)+h_s(1),'--k');
plot(t_plot,y_PI(2,:)+h_s(2),'b'); plot(t_plot,r(2,:)+h_s(2),'-.k')
hold off; grid on; ; title('PI Controller'); ylim([80 140])
subplot(2,3,5); hold on
plot(t_plot,u_PI(1,:)+300,'r'); plot(t_plot,u_PI(2,:)+300,'b');
grid on; ; ylim([0 600])

subplot(2,3,3); hold on
plot(t_plot,y_PID(1,:)+h_s(1),'r'); plot(t_plot,r(1,:)+h_s(1),'--k');
plot(t_plot,y_PID(2,:)+h_s(2),'b'); plot(t_plot,r(2,:)+h_s(2),'-.k')
hold off; grid on; title('PID Controller'); ylim([80 140])
subplot(2,3,6); hold on
plot(t_plot,u_PID(1,:)+300,'r'); plot(t_plot,u_PID(2,:)+300,'b');
grid on; ; ylim([0 600])

%% Functions

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

function [x11,x12,x21,x22] = MarkovParamters(A,B,C,D,N)
    x11 = zeros(1,N+1); x12 = zeros(1,N+1); x21 = zeros(1,N+1); x22 = zeros(1,N+1); 
    Prod = B;
    for k = 1:N+1
        if k == 1
            x11(1,k) = D(1,1);
            x12(1,k) = D(1,2);
            x21(1,k) = D(2,1);
            x22(1,k) = D(2,2);
        else
            H = C*Prod;
            x11(1,k) = H(1,1);
            x12(1,k) = H(1,2);
            x21(1,k) = H(2,1);
            x22(1,k) = H(2,2);
            Prod = A*Prod;
        end % if
    end % k
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

function [x,y,x_hat,y_hat] = KalmanSta(sys,sys_aug,Q_hat,Q_hat_aug,sigma_R,Rd,T_lin,x_s,v_k,d_k,u)
    A = sys.A; B = sys.B(:,1:2); E = sys.B(:,3:4);C = sys.C;
    A_aug = sys_aug.A; B_aug = sys_aug.B(:,1:2); E_aug = sys_aug.B(:,3:4);C_aug = sys_aug.C;

    x = []; x(:,1) = x_s;
    x_hat = []; x_hat(:,1) = [x_s' ; zeros(size(E,2),1)];
    y = []; y_hat = [];

    % See week 6 slide 30
    R = eye(size(C,1),size(C,1)) * sqrt(sigma_R);
    P = idare(A_aug',C_aug',Q_hat_aug,R);
    Re = C_aug*P*C_aug'+R;
    L = P*C_aug'*inv(Re); % Filter can also use DLQE
    
    for i = 1:length(T_lin)-1
        % Measurements
        x(:,i+1) = A*x(:,i) + B*u + E*d_k(:,i);
        y(:,i) = C*x(:,i) + v_k(:,i);
    
        % Estimates and filtering
        y_hat(:,i) = C_aug*x_hat(:,i);
        e(:,i) = y(:,i) - y_hat(:,i);
        x_hat(:,i+1) = A_aug*x_hat(:,i) + B_aug*u +  L*e(:,i);
    end % i
    x = x(:,1:end-1);
    x_hat = x_hat(:,1:end-1);
end % function

function [x,y,x_hat,y_hat] = KalmanDyn(sys,sys_aug,Q_hat,Q_hat_aug,sigma_R,Rd,T_lin,x_s,v_k,d_k,u)
    A = sys.A; B = sys.B(:,1:2); E = sys.B(:,3:4);C = sys.C;
    A_aug = sys_aug.A; B_aug = sys_aug.B(:,1:2); E_aug = sys_aug.B(:,3:4);C_aug = sys_aug.C;

    x = []; x(:,1) = x_s;
    x_hat = []; x_hat(:,1) = [x_s' ; zeros(size(E,2),1)];
    y = []; y_hat = [];

    % See week 6 slide 30
    R = eye(size(C,1),size(C,1)) * sqrt(sigma_R);
    P = idare(A_aug',C_aug',Q_hat_aug,R);
    
    for i = 1:length(T_lin)-1
        % Measurements
        x(:,i+1) = A*x(:,i) + B*u + E*d_k(:,i);
        y(:,i) = C*x(:,i) + v_k(:,i);
    
        % Estimates and filtering
        y_hat(:,i) = C_aug*x_hat(:,i);
        e(:,i) = y(:,i) - y_hat(:,i);
        
        Re = C_aug*P*C_aug'+R;
        L = P*C_aug'*inv(Re);
        
        x_hat(:,i+1) = A_aug*x_hat(:,i) + B_aug*u +  L*e(:,i);
    
        P = A_aug*(P-L*Re*L')*A_aug' + Q_hat_aug;
    end % i
    x = x(:,1:end-1);
    x_hat = x_hat(:,1:end-1);
end % function

function [x_hat,y_hat,P] = oneStepKalman(sys_aug,Q,sigma_R,P,R,Re,L,x_hat,y,u,type)
    A_aug = sys_aug.A; B_aug = sys_aug.B(:,1:2); C_aug = sys_aug.C;
    
    % Static
    if type == 1
        y_hat = C_aug*x_hat;
        e = y - y_hat;
        x_hat = A_aug*x_hat + B_aug*u' +  L*e;
        P = P;
    end

    % Dynamic
    if type == 2 
        y_hat = C_aug*x_hat;
        e = y - y_hat;
        Re = C_aug*P*C_aug'+R;
        L = P*C_aug'*inv(Re);
        x_hat = x_hat + L*e;
        P = P-L*Re*L';

        x_hat = A_aug*x_hat + B_aug*u';
        P = A_aug*P*A_aug'+Q;
    end
end

function  x_out = f_func(x,u,d,p)
    x_out = zeros(size(x,1),1);
    a1 = p(1,1); a2 = p(2,1); a3 = p(3,1); a4 = p(4,1); A1 = p(5,1); A2 = p(6,1); A3 = p(7,1); A4 = p(8,1); gamma1 = p(9,1); gamma2 = p(10,1); g = p(11,1); rho = p(12,1);
    x_out(1,1) = rho*(gamma1*u(1,1) + a3*sqrt(2*g*x(3,1)/(rho*A3)) - a1*sqrt(2*g*x(1,1)/(rho*A1)));
    x_out(2,1) = rho*(gamma2*u(2,1) + a4*sqrt(2*g*x(4,1)/(rho*A4)) - a2*sqrt(2*g*x(2,1)/(rho*A2)));
    x_out(3,1) = rho*((1-gamma2)*u(2,1) + d(1,1) - a3*sqrt(2*g*x(3,1)/(rho*A3)));
    x_out(4,1) = rho*((1-gamma1)*u(1,1) + d(1,1) - a4*sqrt(2*g*x(4,1)/(rho*A4)));
end

function x_out = sigma_func(x,p)
    x_out = zeros(size(x,1),2);
    rho = p(12,1);
    x_out(1,1)=0;
    x_out(2,1)=0;
    x_out(3,1)=rho;
    x_out(4,2)=rho;
end

function [x_k,P_k,x_k1,P_k1,e_k,Re_k,y_k1] = CDEKF(sigma,P_func,p,x_k1,u,d,Pk1,y_k,w_k,t,i,C)
    A1 = p(5,1); rho = p(12,1);
    y_k1 = x_k1(1:2,:) / (rho*A1);
    e_k = y_k1 - y_k;
    Re_k = C*Pk1*C' + w_k;
    K_k = Pk1*C'*inv(Re_k);
    x_k = x_k1 + K_k*e_k;
    P_k = (eye(size(K_k*C,1))-K_k*C) * Pk1 * (eye(size(K_k*C,1))-K_k*C)' + K_k*w_k*K_k';

    [~,P] = ode15s(P_func,[t(i) t(i+1)],[P_k(:,1) ; P_k(:,2) ; P_k(:,3) ; P_k(:,4) ; x_k],[],sigma,[u ; d],p);
    P_k1 = [P(end,1:4)' ; P(end,5:8)' ; P(end,9:12)' ; P(end,13:16)'];
    x_k1 = [P(end,17:20)'];
end

function [P_dot]=P_func(t,P,sigma,in,p)
    a1 = p(1,1); a2 = p(2,1); a3 = p(3,1); a4 = p(4,1); A1 = p(5,1); A2 = p(6,1); A3 = p(7,1); A4 = p(8,1); gamma1 = p(9,1); gamma2 = p(10,1); g = p(11,1); rho = p(12,1);
    u = in(1:2,1); d = in(3:4,1);
        
    x_dot1=-rho*a1*sqrt(2*g*P(17)/(rho*A1))+rho*a3*sqrt(2*g*P(19)/(rho*A3))+rho*gamma1*u(1,1);
    x_dot2=-rho*a2*sqrt(2*g*P(18)/(rho*A2))+rho*a4*sqrt(2*g*P(20)/(rho*A4))+rho*gamma2*u(2,1);
    x_dot3=-rho*a3*sqrt(2*g*P(19)/(rho*A3))+rho*(1-gamma2)*u(2,1)+rho*d(1,1);
    x_dot4=-rho*a4*sqrt(2*g*P(20)/(rho*A4))+rho*(1-gamma1)*u(1,1)+rho*d(1,1);

    A = [-a1/A1*sqrt(g/(2*P(17)/(rho*A1))), 0, a3/A3*sqrt(g/(2*P(19)/(rho*A3))),0;
        0, -a2/A2*sqrt(g/(2*P(18)/(rho*A2))), 0, a4/A4*sqrt(g/(2*P(20)/(rho*A4)));
        0, 0 ,-a3/A3*sqrt(g/(2*P(19)/(rho*A3))),0;
        0, 0, 0, -a4/A4*sqrt(g/(2*P(20)/(rho*A4)))];

    P_1=[P(1:4,1),P(5:8,1),P(9:12,1),P(13:16,1)];
    P_dot1=A*P_1+P_1*(A)'+sigma*(sigma)';
    P_dot=[P_dot1(:,1);P_dot1(:,2);P_dot1(:,3);P_dot1(:,4);[x_dot1;x_dot2;x_dot3;x_dot4]];
end

function [y,y_nl,u,u_nl] = MPC_Sim(sys,sys_aug,MPC_sys,MPC_sys_soft,Q_hat,Q_hat_aug,sigma_R,t_f,t_R,Ts,inputs,N,cons,cons_InOut,cons_Econ,p,type)
    x = cell2mat(inputs(1,1));  x_s = cell2mat(inputs(1,2)); u = cell2mat(inputs(1,3)); 
    R = cell2mat(inputs(1,4)); v_k = cell2mat(inputs(1,5)); d_k = cell2mat(inputs(1,6));
    h_s = ( x_s(:,1:2)/(p(12)*p(5)) )';
    
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
    H = transpose(Gamma)*Q_z*Gamma + Hs;
    
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

function MPC_sys_soft = Soft_MPCDesign(sys,MPC_sys,N,W,soft_cons)
    % Extract data
    W_t1 = cell2mat(soft_cons(1,1)); W_t2 = cell2mat(soft_cons(1,2)); W_s1 = cell2mat(soft_cons(1,3)); W_s2 = cell2mat(soft_cons(1,4));
    Gamma = cell2mat(MPC_sys(1,3)); Lambda = cell2mat(MPC_sys(1,12)); I_0 = cell2mat(MPC_sys(1,13));
    W_z = W(:,1:2); W_u = W(:,3:4); W_du = W(:,5:6); 

    
    % Determine Ht
    W_t2_bar=kron(eye(N),W_t2);
    H_t=W_t2_bar*W_t2_bar;
        
    % Determine gt
    W_t1_bar=kron(eye(N),W_t1);
    g_t=W_t1_bar*ones(size(W_t1_bar,2),1);
    
    % Set point
    W_z_bar=kron(eye(N),W_z);
    H_z=(W_z_bar*Gamma)'*(W_z_bar*Gamma);
    M_z=-(W_z_bar*Gamma)'*W_z_bar;
    
    % Input to reference input
    W_u_bar=kron(eye(N),W_u);
    H_u=W_u_bar'*W_u_bar;
    M_u=-H_u;
    
    % Input variations
    W_du_bar=kron(eye(N),W_du);
    H_du=(W_du_bar*Lambda)'*(W_du_bar*Lambda);
    M_du=-(W_du_bar*Lambda)'*W_du_bar*I_0;
    
    % Determine H total
    H_top=H_z+H_u+H_u;
    
    % lower bound
    W_s2_bar=kron(eye(N),W_s2);   
    W_s1_bar=kron(eye(N),W_s1);
    H_s=W_s2_bar'*W_s2_bar;
    g_s=W_s1_bar*ones(size(W_s1_bar,2),1);

    MPC_sys_soft = {H_t,H_u,H_du,H_top,H_s,g_t,g_s,M_z,M_u,M_du};
end

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

function u_mpc = Uncon_MPC(MPC_sys,x,R,u_prev)
    % Extract data for design
    M_x0 = cell2mat(MPC_sys(1,6)); M_r = cell2mat(MPC_sys(1,7)); M_u1 = cell2mat(MPC_sys(1,9));
    Hs = cell2mat(MPC_sys(1,10)); H = cell2mat(MPC_sys(1,11));
    g = M_x0*x + M_r*R + M_u1*u_prev;
    [u_mpc] = qpsolver(H+Hs,g,[],[],[],[],[],[]);
end % function

function u_mpc = U_con_MPC(MPC_sys,x,R,u_prev,cons);
    % Extract data for design
    M_x0 = cell2mat(MPC_sys(1,6)); M_r = cell2mat(MPC_sys(1,7)); M_u1 = cell2mat(MPC_sys(1,9));
    Hs = cell2mat(MPC_sys(1,10)); H = cell2mat(MPC_sys(1,11)); Lambda = cell2mat(MPC_sys(1,12)); I_0 = cell2mat(MPC_sys(1,13));
    % Extract data for constraints
    lb = cell2mat(cons(1,1)); ub = cell2mat(cons(1,2)); Delta_min = cell2mat(cons(1,3)); Delta_max = cell2mat(cons(1,4));

    ub_Delta = Delta_max+I_0*u_prev;
    lb_Delta = Delta_min+I_0*u_prev;
    g = M_x0*x + M_r*R + M_u1*u_prev;
    [u_mpc] = qpsolver(H+Hs,g,lb,ub,Lambda,lb_Delta,ub_Delta,[]);
end % function

function u_mpc = InOut_MPC(MPC_sys,MPC_sys_soft,x,R,u_prev,cons_InOut,N);
    % Extract data
    H_t = cell2mat(MPC_sys_soft(1,1)); H_top= cell2mat(MPC_sys_soft(1,4)); H_s = cell2mat(MPC_sys_soft(1,5));
    g_t = cell2mat(MPC_sys_soft(1,6));  g_s = cell2mat(MPC_sys_soft(1,7));
    M_z = cell2mat(MPC_sys_soft(1,8)); M_u = cell2mat(MPC_sys_soft(1,9)); M_du = cell2mat(MPC_sys_soft(1,10));

    phi = cell2mat(MPC_sys(1,1)); Gamma = cell2mat(MPC_sys(1,3)); Lambda = cell2mat(MPC_sys(1,12)); I_0 = cell2mat(MPC_sys(1,13)); 
    
    U_bar = cell2mat(cons_InOut(1,1)); U_min = cell2mat(cons_InOut(1,2)); U_max = cell2mat(cons_InOut(1,3)); 
    D_u_min = cell2mat(cons_InOut(1,4)); D_u_max = cell2mat(cons_InOut(1,5));
    Z_max = cell2mat(cons_InOut(1,6)); Z_min = cell2mat(cons_InOut(1,7));

    % Hbar is H matrix quadprog
    H_bar = kron(diag([1 0 0]),H_top)+kron(diag([0 1 0]),H_s)+kron(diag([0 0 1]),H_t);
    b_k = phi*x;
    c_k = R-b_k;
    g = M_z*c_k+M_u*U_bar+M_du*u_prev;
    
    % gbar is g quadprog
    g_bar = [g ; g_s ; g_t];
    
    % Determine A matrix
    A_bar = kron(diag([1 0 0]),Lambda)+kron(diag([0 1 -1]) , eye(2*N))+kron([0 0 0; 1 0 0; 1 0 0] , Gamma);

    % boundary values on states
    u_low = [U_min ; zeros(size(U_min)) ; zeros(size(U_min))];
    
    u_up = [U_max;Inf*ones(size(U_max));inf(1)*ones(size(U_max))];
    
    b_l_bar = [D_u_min+I_0*u_prev;Z_min-b_k;-Inf*ones(size(D_u_min))];
    
    b_u_bar = [D_u_max+I_0*u_prev;inf(1)*ones(size(D_u_max));Z_max-b_k];
    
    u_mpc = qpsolver(H_bar,g_bar,u_low,u_up,A_bar,b_l_bar,b_u_bar,[]);

end

function u_mpc = Econ_MPC(MPC_sys,x,R,u_prev,cons,cons_Econ)
    % Unpack
    phi = cell2mat(MPC_sys(1,1)); Gamma = cell2mat(MPC_sys(1,3)); Lambda = cell2mat(MPC_sys(1,12)); I_0 = cell2mat(MPC_sys(1,13)); 
    lb = cell2mat(cons(1,1)); ub = cell2mat(cons(1,2)); Delta_min = cell2mat(cons(1,3)); Delta_max = cell2mat(cons(1,4));
    g_u = cons_Econ(:,1); g_v = cons_Econ(:,2);
    
    N = length(g_u)/2;
    
    R_bar = R-phi*x;
    
    ub = [ub ; inf(N*2,1)];
    lb = [lb ; zeros(N*2,1)];
    ub_Delta = [Delta_max+I_0*u_prev ; inf(N*2,1)];
    lb_Delta = [Delta_min+I_0*u_prev ; R_bar];

    A = [Lambda , zeros(N*2,N*2) ; Gamma, eye(N*2)];
    
    g = [g_u ; g_v];    
    u_mpc = qpsolver([],g,lb,ub,A,lb_Delta,ub_Delta,[]);
end

function [x_new] = qpsolver(H,g,low_x,up_x,mat,lb,ub,x_init)
    A_quad = [mat ; -mat];
    bound = [ub ; -lb];
    options=optimoptions("quadprog","Display","none");
    [x_new info] = quadprog(H,g,A_quad,bound,[],[],low_x,up_x,x_init,options);
end

function [u e] = P_con(r,y,u_prev,u_min,u_max,K_P)
    e = r-y;
    u = u_prev + K_P*e ;
    if u(1,:)>u_max
        u(1,:)=u_max;
    elseif u(1,:)<u_min
        u(1,:)=u_min;
    end % end if
    if u(2,:)>u_max
        u(2,:)=u_max;
    elseif u(2,:)<u_min
        u(2,:)=u_min;
    end % end if
end

function [u e] = PI_con(r,y,u_prev,u_min,u_max,K_PI,I,Ts,Ti)
    e = r-y;
    u = u_prev + K_PI*e + I;
    if u(1,:)>u_max
        u(1,:)=u_max;
    elseif u(1,:)<u_min
        u(1,:)=u_min;
    else
        I(1,:) = I(1,:) + (K_PI(1,1)*Ts/Ti*e(1,:));
    end % end if
    if u(2,:)>u_max
        u(2,:)=u_max;
    elseif u(2,:)<u_min
        u(2,:)=u_min;
    else
        I(2,:) = I(2,:) + (K_PI(2,2)*Ts/Ti*e(2));
    end % end if
end

function [u e] = PID_con(r,u_prev,y,y_prev,u_min,u_max,K_P,I,K_I,K_D,dt)
    e = r-y;
    P = K_P*e;
    D = -K_D*(y-y_prev)/dt;
    u = u_prev + P + I + D;
    
    if (u(1,:) >= u_max)
        u(1,:) = u_max;
    elseif (u(1,:) <= u_min)
       u(1,:) = u_min;
    else
        I(1,:) = I(1,:)+ K_I(1,1)*e(1,:)*dt;
    end
    
    if (u(2,:) >= u_max)
        u(2,:) = u_max;
    elseif (u(2,:) <= u_min)
        u(2,:) = u_min;
    else
        I(2,:) = I(2,:)+ K_I(2,2)*e(2,:)*dt;
    end
end

function plotFunc(STATES,TIME,n,StateNames,x_label,y_label,leg,name,col)
    subName = string(StateNames);
    name = string(name);
    leg = string(leg); legTemp = cell(size(STATES,2),1); % Used for legend plotting
    if size(STATES,2) > 1;
        for i = 1:size(STATES,2)
            legTemp{i} = leg(i);
        end % i
    end % if
    
    figure
    for i = 1:n
        if n == 2;
            subplot(1,2,i)
            hold on
            for k = 1:size(STATES,2)
                plot(TIME{1,k},STATES{1,k}(:,i),col(k))
                if k > 1; if i == 1; legend(legTemp,'Location','best'); end; end % Ensure corret plotting of legends
            end % k
            hold off
            title(num2str(subName(i))); xlabel(x_label); ylabel(y_label); sgtitle(name); grid on
        else % else n NOT EQUAL 2
            plot_order = [3 4 1 2]; % Used to shift the order of subplots
            subplot(2,2,plot_order(i))
            hold on
            for k = 1:size(STATES,2)
                   plot(TIME{1,k},STATES{1,k}(:,i),col(k))
               if k > 1; if i == 3; legend(legTemp,'Location','best'); end; end % Ensure corret plotting of legends
            end % k
            hold off
            title(num2str(subName(plot_order(i)))); xlabel(x_label); ylabel(y_label); sgtitle(name); grid on
        end% n
    end % i
end % Function
function [W,t,dW] = ScalarStdWienerProcess(T,N,Ns,seed)
% ScalarStdWienerProcess Ns realizations of a scalar std Wiener process
%
% Syntax: [W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed)
% W : Standard Wiener process in [0,T]
% Tw : Time points
% dW : White noise used to generate the Wiener process
%
% T : Final time
% N : Number of intervals
% Ns : Number of realizations
% seed : To set the random number generator (optional)

if nargin == 4 
rng(seed);
end
%nargin er en indbygget funktion i MATLAB, der returnerer antallet af 
% inputargumenter, der er givet til den aktuelle funktion.

% hver af disse linjer forklares i Scalar_Standard_Brownian_Motion_Standard_Wiener
dt = T/N;
dW = [zeros(1,Ns); sqrt(dt)*randn(N,Ns)];
W = [cumsum(dW,1)];
t = (0:dt:T)';
end


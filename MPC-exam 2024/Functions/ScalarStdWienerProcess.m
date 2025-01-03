function [W, t, dW] = ScalarStdWienerProcess(T, N, Ns, seed)
% Generates Ns realizations of a scalar standard Wiener process.
% Set random number generator seed if specified
if nargin == 4
    rng(seed);  % Initialize random number generator for reproducibility
end
% Compute time step size
dt = T / N;    
% Generate white noise increments
dW = [zeros(1, Ns); sqrt(dt) * randn(N, Ns)];  
% Compute the Wiener process using cumulative sum
W = cumsum(dW, 1);  
% Create a time vector from 0 to T with step size dt
t = (0:dt:T)';  
end

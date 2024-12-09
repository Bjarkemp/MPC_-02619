function [A, B, C, G, Gw] = linearized_models(xs, At, at, rho, gamma1, gamma2, g, noise_type, d0, sigma, R)
    % LINEARIZED_MODELS computes linearized continuous-time models for different noise types.
    %
    % Inputs:
    %   xs        - Steady-state masses in the tanks.
    %   At        - Tank cross-sectional areas.
    %   at        - Outlet cross-sectional areas.
    %   rho       - Liquid density.
    %   gamma1    - Split fraction for Tank 1.
    %   gamma2    - Split fraction for Tank 2.
    %   g         - Gravitational constant.
    %   noise_type- 'stochastic' or 'brownian'.
    %   d0        - Base disturbance vector for Brownian motion.
    %   sigma     - Covariance matrix for disturbances (Brownian motion).
    %   R         - Covariance matrix for measurement noise.
    %
    % Outputs:
    %   A, B, C   - State-space matrices.
    %   G         - Disturbance propagation matrix.
    %   Gw        - Measurement noise matrix.

    % Compute state-space matrices for the base system
    [A, B, C, ~, G, Gw] = system_matrices(xs, At, at, rho, gamma1, gamma2, g);

    % Adjust for noise type
    if strcmp(noise_type, 'stochastic')
        % Stochastic disturbances
        d = [randn(1, length(xs)); randn(1, length(xs))] * 2;  % Gaussian noise
        G = G * d;  % Adjust disturbance propagation

    elseif strcmp(noise_type, 'brownian')
        % Brownian motion disturbances
        G = G * sigma;  % Scale disturbance propagation by covariance
    end

    % Augment measurement noise
    Gw = eye(size(A, 1));  % Independent noise for each state
end

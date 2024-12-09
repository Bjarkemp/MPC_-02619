function [A, B, C, Cz, G, Gw] = system_matrices(xs, At, at, rho, gamma1, gamma2, g)
    % SYSTEM_MATRICES computes the state-space matrices for a tank system.
    %
    % Inputs:
    %   xs     - Vector of masses in the tanks (steady-state values).
    %   At     - Cross-sectional area of the tanks.
    %   at     - Cross-sectional area of the outlets.
    %   rho    - Density of the liquid.
    %   gamma1 - Split fraction for input to Tank 1.
    %   gamma2 - Split fraction for input to Tank 2.
    %   g      - Gravitational acceleration constant.
    %
    % Outputs:
    %   A      - State matrix.
    %   B      - Input matrix.
    %   C      - Output matrix for all tanks.
    %   Cz     - Output matrix for first two tanks.
    %   G      - Disturbance matrix.
    %   Gw     - Noise matrix.

    % Convert mass to height
    hs = mass_to_height(xs, At, rho);

    % Compute time constants
    Tl = (At ./ at) .* sqrt(2 * hs / g);

    % Define the A matrix (state dynamics)
    A = [-1/Tl(1),  0,        1/Tl(3),  0;
          0,      -1/Tl(2),  0,        1/Tl(4);
          0,       0,       -1/Tl(3),  0;
          0,       0,        0,       -1/Tl(4)];

    % Define the B matrix (input dynamics)
    B = [rho * gamma1,        0;
         0,             rho * gamma2;
         0,             rho * (1 - gamma2);
         rho * (1 - gamma1),  0];

    % Define the C matrix (output dynamics for all tanks)
    C = diag(1 ./ (rho * At));

    % Output matrix for first two tanks
    Cz = C(1:2, :);

    % Disturbance matrix
    G = [0,    0;
         0,    0;
         rho,  0;
         0,  rho];

    % Noise matrix (identity for independent noise in each state)
    Gw = eye(4);
end


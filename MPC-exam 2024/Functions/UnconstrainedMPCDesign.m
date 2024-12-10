function MPC_sys = UnconstrainedMPCDesign(A, B, C, Q, S, N)
    % UNCONSTRAINEDMPCDESIGN designs an Unconstrained MPC for a discrete-time system.
    %
    % Inputs:
    %   A, B, and C - Discrete-time state-space system (struct with A, B, C matrices)
    %   Q   - State tracking weight matrix
    %   S   - Input weight matrix
    %   N   - Prediction horizon
    %
    % Output:
    %   MPC_sys - Struct containing MPC matrices (phi, Gamma, H, g, etc.)

    % Generate prediction matrices (phi and Gamma)
    phi = generate_phi(A, C, N);      % State-to-output mapping
    Gamma = generate_Gamma(A, B, C, N); % Input-to-output mapping

    % Define Q_z (weight on predicted outputs)
    Q_z = kron(eye(N), Q);  % Block diagonal weight matrix

    % Define H (quadratic term in cost function)
    H = Gamma' * Q_z * Gamma + kron(eye(N), S); % Prediction weight + control weight

    % Define M_x0 (linear term related to initial state)
    M_x0 = Gamma' * Q_z * phi;

    % Define M_r (linear term related to reference tracking)
    M_r = -Gamma' * Q_z;

    % Return MPC system matrices
    MPC_sys = struct('phi', phi, 'Gamma', Gamma, 'H', H, 'M_x0', M_x0, 'M_r', M_r);
end

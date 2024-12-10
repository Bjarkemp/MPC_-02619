function MPC_sys = UnconstrainedMPCDesign(A, B, C, Q, S, Ph)
    % UNCONSTRAINEDMPCDESIGN designs an Unconstrained MPC for a discrete-time system.
    %
    % Inputs:
    %   A, B, and C - Discrete-time state-space system (struct with A, B, C matrices)
    %   Q   - State tracking weight matrix
    %   S   - Input weight matrix
    %   Ph   - Prediction horizon
    %
    % Output:
    %   MPC_sys - Struct containing MPC matrices (phi, Gamma, H, g, etc.)

    % Generate prediction matrices (phi and Gamma)
    phi = generate_phi(A, C, Ph);      % State-to-output mapping
    Gamma = generate_Gamma(A, B, C, Ph); % Input-to-output mapping

    % Define Q_z (weight on predicted outputs)
    Q_z = kron(eye(Ph), Q);  % Block diagonal weight matrix

    % Generate Hs
    Hs = zeros(Ph*2);
    % Create arrow for matrix dimensions
    pil_Hs_row = [1:size(S,1)]; pil_Hs_col = pil_Hs_row;
    for i = 1:Ph
        while (pil_Hs_row(end) < Ph*2)
            Hs(pil_Hs_row,pil_Hs_col) = S*2; % Set S in diagonal
            Hs(pil_Hs_row,pil_Hs_col+size(S,1)) = -S;
            Hs(pil_Hs_row+size(S,1),pil_Hs_col) = -S;
            % Update arrows
            pil_Hs_row = pil_Hs_row + size(S,1); pil_Hs_col = pil_Hs_row; 
        end % while
        if i == Ph
            Hs(pil_Hs_row,pil_Hs_col) = S;
        end % if i == Ph
        % Update arrows
        pil_Hs_row = [1:size(S,1)] + size(S,1)*i; pil_Hs_col = pil_Hs_row;
    end % i

    % Define H (quadratic term in cost function)
    H = Gamma' * Q_z * Gamma + Hs; % Prediction weight + control weight

    % Define M_x0 (linear term related to initial state)
    M_x0 = Gamma' * Q_z * phi;

    % Define M_r (linear term related to reference tracking)
    M_r = -Gamma' * Q_z;

    % Return MPC system matrices
    MPC_sys = struct('phi', phi, 'Gamma', Gamma, 'H', H, 'M_x0', M_x0, 'M_r', M_r);
end

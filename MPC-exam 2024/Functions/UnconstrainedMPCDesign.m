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

% Initialize Hs matrix
Hs = zeros(2 * Ph);

% Initialize row and column indices
row_idx = 1:size(S, 1);
col_idx = row_idx;

% Loop through the phases
for phase = 1:Ph
    % Fill in the submatrices for the current phase
    while row_idx(end) < 2 * Ph
        % Set values in the main diagonal and submatrices
        Hs(row_idx, col_idx) = 2 * S;
        Hs(row_idx, col_idx + size(S, 1)) = -S;
        Hs(row_idx + size(S, 1), col_idx) = -S;

        % Update indices for the next set of blocks
        row_idx = row_idx + size(S, 1);
        col_idx = row_idx;
    end

    % Handle the last phase separately
    if phase == Ph
        Hs(row_idx, col_idx) = S;
    end

    % Reset row and column indices for the next phase
    row_idx = (1:size(S, 1)) + size(S, 1) * phase;
    col_idx = row_idx;
end

    % Define H (quadratic term in cost function)
    H = Gamma' * Q_z * Gamma + Hs; % Prediction weight + control weight

    % Define M_x0 (linear term related to initial state)
    M_x0 = Gamma' * Q_z * phi;

    % Define M_r (linear term related to reference tracking)
    M_r = -Gamma' * Q_z;
    
    % Define M_u1 (linear term related to control constraints)
    M_u1 = zeros(size(Hs,1),size(S,1)); 
    M_u1(1:size(S,1),1:size(S,1)) = -S;

    % Determine the input bounds
    Lambda = eye(Ph); 
    Lambda = kron(eye(2), Lambda); % Expand Gamma for block structure

    % Update the Lambda matrix to include bounds
    for idx = 1:2*Ph
        Lambda(idx+2:idx+3, idx:idx+1) = -eye(2);
    end
    
    % Trim the last rows and columns to finalize the matrix
    Lambda = Lambda(1:end-3, 1:end-1);
    
    % Initialize the I_base matrix with zeros
    I0 = zeros(Ph * size(S, 1), size(S, 2));
    
    % Set the initial identity block
    I0(1:size(S, 1), 1:size(S, 2)) = eye(size(S, 1));

    % Return MPC system matrices
    MPC_sys = struct('phi', phi, 'Gamma', Gamma, 'H', H, 'Hs', Hs, 'M_x0', M_x0, 'M_r', M_r, 'M_u1', M_u1, 'Lambda', Lambda, 'I0', I0);
end

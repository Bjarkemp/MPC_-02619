function phi_w = generate_phi_w(A, C, E, N)
    % GENERATE_PHI_W Computes the phi_w matrix for MPC.
    %
    % Inputs:
    %   A - State matrix
    %   C - Output matrix
    %   E - Disturbance matrix
    %   N - Prediction horizon
    %
    % Output:
    %   phi_w - Mapping from disturbances to predicted outputs

    phi_w = zeros(N * size(C, 1), size(E, 2)); % Initialize matrix
    pil_phi_w = 1:size(C, 1); % Row index tracker

    for i = 1:N
        phi_w(pil_phi_w, :) = C * A^i * E; % Compute C * A^i * E
        pil_phi_w = pil_phi_w + size(C, 1); % Update row indices
    end
end

function phi = generate_phi(A, C, N)
    % GENERATE_PHI Computes the phi matrix for MPC.
    %
    % Inputs:
    %   A - State matrix
    %   C - Output matrix
    %   N - Prediction horizon
    %
    % Output:
    %   phi - Mapping from initial state to predicted outputs
    
    phi = zeros(N * size(C, 1), size(A, 2)); % Initialize matrix
    pil_phi = 1:size(C, 1); % Row index tracker

    for i = 1:N
        phi(pil_phi, :) = C * A^i; % Compute C * A^i and place it
        pil_phi = pil_phi + size(C, 1); % Update row indices
    end
end

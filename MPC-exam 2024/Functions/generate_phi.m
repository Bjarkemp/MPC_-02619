function phi = generate_phi(A, C, N)
    % GENERATE_PHI calculates the phi matrix for the MPC horizon
    %
    % Inputs:
    %   A - System state matrix
    %   C - System output matrix
    %   N - Prediction horizon
    %
    % Outputs:
    %   phi - Matrix mapping initial state to predicted outputs over the horizon

    % Initialize phi
    phi = zeros(N * size(C, 1), size(A, 2));

    % Create phi for N iterations
    for i = 1:N
        % Compute phi row block for iteration i
        phi((i-1)*size(C, 1) + 1:i*size(C, 1), :) = C * A^i;
    end
end

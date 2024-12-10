function phi = generate_phi(A, C, Ph)
    % GENERATE_PHI calculates the phi matrix for the MPC horizon
    %
    % Inputs:
    %   A - System state matrix
    %   C - System output matrix
    %   Ph - Prediction horizon
    %
    % Outputs:
    %   phi - Matrix mapping initial state to predicted outputs over the horizon

    % Initialize phi
    phi = zeros(Ph * size(C, 1), size(A, 2));

    % Create phi for Ph iterations
    for i = 1:Ph
        % Compute phi row block for iteration i
        phi((i-1)*size(C, 1) + 1:i*size(C, 1), :) = C * A^i;
    end
end

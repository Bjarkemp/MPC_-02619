function Gamma = generate_Gamma(A, B, C, N)
    % GENERATE_GAMMA calculates the Gamma matrix for the MPC horizon
    %
    % Inputs:
    %   A - System state matrix
    %   B - System input matrix
    %   C - System output matrix
    %   N - Prediction horizon
    %
    % Outputs:
    %   Gamma - Matrix mapping input sequence to predicted outputs over the horizon

    % Initialize Gamma
    Gamma = zeros(N * size(C, 1), N * size(B, 2));

    % Populate Gamma for N iterations
    for i = 1:N
        for j = 1:i
            % Compute Gamma block at (i, j)
            row_start = (i-1)*size(C, 1) + 1;
            row_end = i*size(C, 1);
            col_start = (j-1)*size(B, 2) + 1;
            col_end = j*size(B, 2);

            Gamma(row_start:row_end, col_start:col_end) = C * A^(i-j) * B;
        end
    end
end

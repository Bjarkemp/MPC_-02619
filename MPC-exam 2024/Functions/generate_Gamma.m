function Gamma = generate_Gamma(A, B, C, Ph)
    % GENERATE_GAMMA calculates the Gamma matrix for the MPC horizon
    %
    % Inputs:
    %   A - System state matrix
    %   B - System input matrix
    %   C - System output matrix
    %   Ph - Prediction horizon
    %
    % Outputs:
    %   Gamma - Matrix mapping input sequence to predicted outputs over the horizon

    % Initialize Gamma
    Gamma = zeros(Ph * size(C, 1), Ph * size(B, 2));

    % Populate Gamma for Ph iterations
    for i = 1:Ph
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

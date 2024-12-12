function Gamma = generate_Gamma_d(A, E, C, N)
    % GENERATE_GAMMA Computes the Gamma matrix for MPC as shown in the image.
    %
    % Inputs:
    %   A - State matrix
    %   E - Disturbance matrix
    %   C - Output matrix
    %   N - Prediction horizon
    %
    % Output:
    %   Gamma - Block-lower triangular mapping from inputs to predicted outputs

    % Dimensions
    n_outputs = size(C, 1);  % Number of outputs
    n_inputs = size(E, 2);   % Number of inputs

    % Initialize Gamma
    Gamma = zeros(N * n_outputs, N * n_inputs);

    % Fill the block-lower triangular structure
    for row_block = 1:N
        for col_block = 1:row_block
            % Compute H_i = C * A^(row_block - col_block) * B
            block = C * A^(row_block - col_block) * E;
            
            % Row and column indices for placing the block
            row_indices = (row_block-1) * n_outputs + (1:n_outputs);
            col_indices = (col_block-1) * n_inputs + (1:n_inputs);
            
            % Place the block in the appropriate location in Gamma
            Gamma(row_indices, col_indices) = block;
        end
    end
end

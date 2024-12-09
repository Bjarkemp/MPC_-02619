function G_tf = compute_transfer_functions(A, B, C, D)
    % COMPUTE_TRANSFER_FUNCTIONS calculates the transfer functions for a linearized model
    %
    % Inputs:
    %   A, B, C, D - State-space matrices
    %
    % Outputs:
    %   G_tf - Cell array of transfer functions (one per input-output pair)
    
    % Create state-space model
    sys_ss = ss(A, B, C, D);  
    
    % Convert state-space to transfer function
    G_tf = tf(sys_ss);  % G_tf is a matrix of transfer functions
end
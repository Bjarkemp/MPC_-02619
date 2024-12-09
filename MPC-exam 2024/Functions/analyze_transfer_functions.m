function [gains, time_constants] = analyze_transfer_functions(G_tf)
    % ANALYZE_TRANSFER_FUNCTIONS calculates steady-state gains and time constants
    % for a transfer function matrix.
    %
    % Inputs:
    %   G_tf - Transfer function matrix
    %
    % Outputs:
    %   gains - Steady-state gains for each input-output pair
    %   time_constants - Time constants for each input-output pair

    [num, den] = tfdata(G_tf, 'v');  % Get numerator and denominator arrays
    n_outputs = size(G_tf, 1);
    n_inputs = size(G_tf, 2);
    
    % Preallocate matrices
    gains = zeros(n_outputs, n_inputs);  
    time_constants = zeros(n_outputs, n_inputs);
    
    for i = 1:n_outputs
        for j = 1:n_inputs
            % Check if the transfer function is zero
            if all(num{i, j} == 0)
                % Zero transfer function, set gain and time constant to NaN
                gains(i, j) = NaN;
                time_constants(i, j) = NaN;
                continue;
            end
            
            % Compute steady-state gain: DC gain = numerator constant / denominator constant
            gains(i, j) = num{i, j}(end) / den{i, j}(end);
            
            % Compute dominant time constant
            if length(den{i, j}) > 2
                % If second-order or higher, calculate approximate dominant time constant
                roots_den = roots(den{i, j});  % Find poles of the system
                dominant_pole = min(abs(roots_den));  % Smallest pole gives largest time constant
                time_constants(i, j) = 1 / dominant_pole;  % Time constant is inverse of pole
            else
                % For first-order systems: tau = coefficient of s / constant term
                time_constants(i, j) = den{i, j}(end-1) / den{i, j}(end);
            end
        end
    end
end

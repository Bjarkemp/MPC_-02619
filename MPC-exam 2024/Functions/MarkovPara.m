function [m11, m12, m21, m22] = ComputeMarkov(A_mat, B_mat, C_mat, D_mat, steps)
    % Initialize output matrices
    m11 = zeros(1, steps + 1); 
    m12 = zeros(1, steps + 1); 
    m21 = zeros(1, steps + 1); 
    m22 = zeros(1, steps + 1); 
    
    % Initialize product term
    current_prod = B_mat;
    
    for idx = 1:(steps + 1)
        if idx == 1
            % Assign initial values from matrix D
            m11(1, idx) = D_mat(1, 1);
            m12(1, idx) = D_mat(1, 2);
            m21(1, idx) = D_mat(2, 1);
            m22(1, idx) = D_mat(2, 2);
        else
            % Calculate and update output matrices
            result = C_mat * current_prod;
            m11(1, idx) = result(1, 1);
            m12(1, idx) = result(1, 2);
            m21(1, idx) = result(2, 1);
            m22(1, idx) = result(2, 2);
            
            % Update product term
            current_prod = A_mat * current_prod;
        end
    end
end

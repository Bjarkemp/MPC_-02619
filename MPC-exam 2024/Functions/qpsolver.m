function solution = qpsolver(H, g, lower_bound, upper_bound, constraint_mat, constraint_lb, constraint_ub, initial_guess)
    % SOLVE_QP Solves a quadratic programming problem using quadprog.
    %
    % Inputs:
    %   H - Hessian matrix for the quadratic term (n×n)
    %   f - Linear term vector (n×1)
    %   lower_bound - Lower bound on the solution variables
    %   upper_bound - Upper bound on the solution variables
    %   constraint_mat - Matrix for linear inequality constraints
    %   constraint_lb - Lower bound for linear inequality constraints
    %   constraint_ub - Upper bound for linear inequality constraints
    %   initial_guess - Initial guess for the optimization
    %
    % Output:
    %   solution - Solution vector to the QP problem

    % Combine inequality constraints into a single matrix and vector
    combined_constraints = [constraint_mat; -constraint_mat];
    combined_bounds = [constraint_ub; -constraint_lb];
    
    % Set up optimization options
    qp_options = optimoptions("quadprog", "Display", "none");
    
    % Solve the quadratic programming problem
    [solution, ~] = quadprog(H, g, combined_constraints, combined_bounds, [], [], lower_bound, upper_bound, initial_guess, qp_options);
end



% function [x_new] = qpsolver(H,g,low_x,up_x,mat,lb,ub,x_init)
%     A_quad = [mat ; -mat];
%     bound = [ub ; -lb];
%     options=optimoptions("quadprog","Display","none");
%     [x_new info] = quadprog(H,g,A_quad,bound,[],[],low_x,up_x,x_init,options);
% end
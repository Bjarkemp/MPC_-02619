function [x, info] = qpsolver(H, g, l, u, A, bl, bu, xinit)
    % QP Solver Interface
    % Solves the convex quadratic program:
    %   minimize (1/2)*x'*H*x + g'*x
    %   subject to:
    %       l <= x <= u
    %       bl <= A*x <= bu
    %
    % Inputs:
    %   H     - Quadratic cost matrix
    %   g     - Linear cost vector
    %   l     - Lower bounds on decision variables
    %   u     - Upper bounds on decision variables
    %   A     - Constraint matrix
    %   bl    - Lower bounds on linear constraints
    %   bu    - Upper bounds on linear constraints
    %   xinit - Initial guess for the solution (not used in quadprog directly)
    %
    % Outputs:
    %   x     - Optimal solution vector
    %   info  - Solver information (e.g., exit flag, iterations, etc.)

    % Ensure inputs are consistent with quadprog requirements
    % Convert bounds to the inequality form required by quadprog:
    % A_eq*x = b_eq (no equality constraints in this case)
    % A_ineq*x <= b_ineq (combining variable bounds and linear constraints)

    % Combine variable bounds into inequality constraints
    % l <= x <= u becomes:
    %   [I; -I] * x <= [u; -l]
    n = length(g);  % Number of variables
    A_ineq = [eye(n); -eye(n)];
    b_ineq = [u; -l];

    % Add linear constraints: bl <= A*x <= bu
    if ~isempty(A)
        A_ineq = [A_ineq; A; -A];
        b_ineq = [b_ineq; bu; -bl];
    end

    % Call quadprog to solve the QP
    options = optimoptions('quadprog', 'Display', 'off'); % Suppress output
    [x, fval, exitflag, output] = quadprog(H, g, A_ineq, b_ineq, [], [], [], [], xinit, options);

    % Pack solver information
    info.exitflag = exitflag; % Exit condition of the solver
    info.fval = fval;        % Final value of the objective function
    info.output = output;    % Detailed solver output
end


% Example data (Adjusted for consistency)
H = eye(8);                      % Quadratic cost (identity for simplicity)
g = [-2; -5; 0; 0; 0; 0; 0; 0];  % Linear cost (adjusted for dimension)
l = zeros(8, 1);                 % Lower bounds (all variables >= 0)
u = ones(8, 1);                  % Upper bounds (all variables <= 1)

% Constraint matrix and bounds
M = [
    1.1108,    0,  -0.2389,    0,  -0.0699,    0,  -1.0632,    0;
         0, 1.0843,       0, -0.1611,       0, -0.0346,       0, -0.7408;
         0,      0,   1.2270,      0,   0.9914,      0,  10.0699,      0;
         0,      0,        0, 1.1548,        0, 0.7019,        0, 10.0346;
         0,      0,        0,      0,   0.9003,      0,        0,      0;
         0,      0,        0,      0,        0, 0.9222,        0,      0;
         0,      0,        0,      0,   0.1753,      0,   0.8150,      0;
         0,      0,        0,      0,        0, 0.1287,        0,  0.8659;
];
                                  % Using M as the constraint matrix
bl = -ones(8, 1);                % Lower constraint bounds (adjust as needed)
bu = ones(8, 1);                 % Upper constraint bounds (adjust as needed)

% Initial guess
xinit = 0.5 * ones(8, 1);        % Initial guess (all variables start at 0.5)

% Call the solver
[x, info] = qpsolver(H, g, l, u, M, bl, bu, xinit);

% Display results
disp('Optimal solution:');
disp(x);
disp('Solver info:');
disp(info);

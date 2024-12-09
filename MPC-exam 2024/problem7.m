% Example data
H = [2, 0; 0, 2];          % Quadratic cost
g = [-2; -5];              % Linear cost
l = [0; 0];                % Lower bounds
u = [1; 1];                % Upper bounds
A = [1, 2; 2, 1];          % Constraint matrix
bl = [1; 1];               % Lower constraint bounds
bu = [2; 2];               % Upper constraint bounds
xinit = [];                % Initial guess

% Call the solver
[x, info] = qpsolver(H, g, l, u, A, bl, bu, xinit);

% Display results
disp('Optimal solution:');
disp(x);
disp('Solver info:');
disp(info);
% Load the data from a CSV file
data = csvread('databp.csv');  % Adjust as needed

% Initial guess
initial_guess = [0.4, 0.2, 0.25, 0.2, 0.3, 0.35, 0.15, 0.2, 0.4, 0.3, 0.2];

% Define linear inequality constraints
% A_ineq has 4 rows (one for each constraint) and 11 columns (for each parameter)
A_ineq = zeros(4, 11);

% Constraint for first row of A
A_ineq(1, [2, 3]) = 1;  % alpha_11 + alpha_12 < 1

% Constraint for second row of A
A_ineq(2, [4, 5]) = 1;  % beta_11 + beta_12 < 1

% Constraint for first row of B
A_ineq(3, [7, 8]) = 1;  % alpha_21 + alpha_22 < 1

% Constraint for second row of B
A_ineq(4, [9, 10]) = 1;  % beta_21 + beta_22 < 1

% b_ineq is the right-hand side of the inequalities (less than 1)
b_ineq = [1; 1; 1; 1];

% Bounds for each parameter (lower and upper bounds)
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1];  % Lower bounds
ub = [2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1];   % Upper bounds

% Run optimization with linear constraints
[theta_opt, fval] = fmincon(@(theta) neg_ll(theta, data), initial_guess, A_ineq, b_ineq, [], [], lb, ub, []);  % No equality constraints or nonlinear constraints
disp('Optimized parameters:');
disp(theta_opt);
disp('Minimum log-likelihood:');
disp(fval);
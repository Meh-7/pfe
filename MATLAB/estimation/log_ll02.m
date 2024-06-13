% Load the data from a CSV file
data = csvread('databp2.csv');  % Adjust as needed

% Initial guess
initial_guess = [0.9, 0.35, 0.3, 0.4, 0.35, 0.35, 0.4];


% Bounds for each parameter (lower and upper bounds)
lb = [0, 0, 0, 0, 0, 0, -1];  % Lower bounds
ub = [1.5, 1, 1, 1.5, 1, 1, 1];   % Upper bounds

% Nonlinear constraint function 
nonlincon = @(theta) constraint_delta(theta, data);
% Run optimization with linear and nonlinear constraints
[theta_opt, fval] = fmincon(@(theta) neg_ll_diag(theta, data), initial_guess, [], [], [], [], lb, ub, nonlincon); 

disp('Optimized parameters:');
disp(theta_opt);
disp('Minimum log-likelihood:');
disp(fval);
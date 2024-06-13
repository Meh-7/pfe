% Load the data from CSV
data = csvread('databp.csv');  % Adjust the file name and location as needed
% Initial guess for optimization based on the new structure of theta
initial_guess = [0.4, 0.2, 0.25, 0.2, 0.3, 0.35, 0.15, 0.2, 0.4, 0.3, 0.2];

% Run the optimization with fmincon
lb=[0,0,0,0,0,0,0,0,0,0,0]; 
ub=[1,1,1,1,1,1,1,1,1,1,1];
theta_opt = fmincon(@(theta) neg_ll(theta, data), initial_guess,[], [], [], [], lb, ub, []);

% Display results
disp('Optimized parameters:');
disp(theta_opt);
disp('Minimum log-likelihood:');
disp(fval)
% Load the data from a CSV file
data = csvread('databp.csv'); 
% initial guess for the solution
theta0 = [0.4, 0.2, 0.25, 0.2, 0.3, 0.35, 0.15, 0.2, 0.4, 0.3, 0.2];

% solve the problem
[solution,fval,exitflag,output] = minimize_cost_function(theta0)
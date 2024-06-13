% Load the data from the CSV file
data = table2array(readtable('databp2.csv', 'HeaderLines', 0)); 

% Set initial parameter values
theta0 = [0.8; 0.3; 0.3; 1; 0.3; 0.3; 0.3]; 

% Define lower and upper bounds for the parameters

lb = [0; 0; 0; 0; 0; 0; -1];  % All parameters except delta have lower bound 0
ub = [1; 1; 1; 200; 1; 1; 1]; % Upper bounds, with 100 for omega1 and omega2

% Set fmincon options
options = optimoptions('fmincon', 'Display', 'iter');

% Estimate the parameters using fmincon
theta_hat = fmincon(@(theta) bivariateINGARCHnll(theta, data), theta0, [], [], [], [], lb, ub, [], options);

% Display the estimated parameters
disp('Estimated Parameters:');
disp(theta_hat);
% Load the data from the CSV file
data = table2array(readtable('databp2.csv'));

% Define the objective function for fmincon (this is your provided function)

% Number of replications
num_replications = 200;

% Store estimated parameters from each replication
theta_hat_all = zeros(num_replications, 7);
%theta0 = rand(7, num_replications);
% Loop for multiple replications
for i = 1:num_replications
    % Randomly select starting values for theta0 from a uniform distribution
    theta0 = rand(7, 1);  % Values between 0 and 1
    %theta0(1) = 0 + (2-0).*rand(1,1);
    %theta0(4) = 0 + (2-0).*rand(1,1);
    %theta0(7) = -2 + (2+2).*rand(1,1);
    
    % Define lower and upper bounds for the parameters
    lb = [0; 0; 0; 0; 0; 0; -1];  
    ub = [2; 1; 1; 2; 1; 1; 1]; 

    % Set fmincon options
    options = optimoptions('fmincon', 'Display', 'off'); % Suppress output
    
    % Estimate the parameters using fmincon
    theta_hat = fmincon(@(theta) bivariateINGARCHnll(theta, data), theta0, [], [], [], [], lb, ub, [], options);

    % Store the estimated parameters
    theta_hat_all(i, :) = theta_hat';
end

% Calculate the mean and standard deviation for each parameter
mean_theta_hat = mean(theta_hat_all);
std_theta_hat = std(theta_hat_all);

% Display results
disp('Mean of Estimated Parameters:');
disp(mean_theta_hat);
disp('Standard Deviation of Estimated Parameters:');
disp(std_theta_hat);
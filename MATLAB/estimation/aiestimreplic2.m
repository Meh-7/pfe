% Load the data from the CSV file
data = table2array(readtable('databp2.csv'));

% Define the objective function for fmincon (this is your provided function)
objective = @(theta) bivariateINGARCHnll(theta, data);

% Number of replications
num_replications = 200;

% Store estimated parameters from each replication
theta_hat_all = zeros(num_replications, 7);

% Loop for multiple replications
for i = 1:num_replications
    % Randomly select starting values for theta0 from a wider uniform distribution
    theta0 = [2*rand(1, 1); rand(2, 1); 2*rand(1, 1); rand(2, 1); -2 + 4*rand(1, 1)];
    
    % Debugging: display initial values
    fprintf('Initial theta0 values for replication %d: %s\n', i, mat2str(theta0));
    
    % Define lower and upper bounds for the parameters
    lb = [0; 0; 0; 0; 0; 0; -1];  
    ub = [2; 1; 1; 2; 1; 1; 1]; 

    % Set fmincon options
    options = optimoptions('fmincon', 'Display', 'off', 'MaxFunEvals', 10000, 'MaxIter', 1000, 'TolFun', 1e-6); % Adjusted settings
    
    try
        % Estimate the parameters using fmincon
        theta_hat = fmincon(objective, theta0, [], [], [], [], lb, ub, [], options);

        % Store the estimated parameters
        theta_hat_all(i, :) = theta_hat';
    catch ME
        % Debugging: display error message
        fprintf('Error in replication %d: %s\n', i, ME.message);
    end
end

% Calculate the mean and standard deviation for each parameter
mean_theta_hat = mean(theta_hat_all);
std_theta_hat = std(theta_hat_all);

% Display results
disp('Mean of Estimated Parameters:');
disp(mean_theta_hat);
disp('Standard Deviation of Estimated Parameters:');
disp(std_theta_hat);
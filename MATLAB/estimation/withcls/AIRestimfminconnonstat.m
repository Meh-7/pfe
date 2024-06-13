% Load the data from the CSV file
data = table2array(readtable('databp_M3_500a.csv'));

% Number of replications
num_replications = 200;

% Store estimated parameters from each replication
theta_hat = zeros(num_replications, 11); % 11 parameters now
theta0 = zeros(num_replications, 11);
% Define lower and upper bounds for the parameters (updated for 11 parameters)
    lb = [0.1; 0; 0; 0; 0.1; 0; 0; 0; 0; 0; -2];  % Lower bounds
    ub = [2; 1; 1; 1; 1; 2; 1; 1; 1; 1; 2];  % Upper bounds
    
% Loop for multiple replications
for i = 1:num_replications
    theta0(i,:) = lb + rand(size(lb)).* (ub - lb);
    % Randomly perturb the ranges for generating initial theta0
    %lb_perturbed = lb - rand(size(lb)) * 0.1;
    %ub_perturbed = ub + rand(size(ub)) * 0.1;
    
    % Randomly select starting values for theta0 
    %theta0(i,:) = lb_perturbed + rand(size(lb_perturbed)) .* (ub_perturbed - lb_perturbed);

    % Ensure theta0 remains within bounds
    %theta0(i,:) = max(lb.', theta0(i,:));
    %theta0(i,:) = min(ub.', theta0(i,:));

    % Set fmincon options
    options = optimset('fmincon');
    options.Display = 'off';
    options.GradObj = 'on';

    % Estimate the parameters using fmincon
    [theta_hat(i,:), fval, exitflag, output] = fmincon(@(theta) RbivariateINGARCHnll2(theta, data,3), ...
                                                     theta0(i,:), [], [], [], [], lb, ub, [], options);
end

% Calculate mean and standard deviation for each parameter
mean_theta_hat = mean(theta_hat);
std_theta_hat = std(theta_hat);

% Display results
disp('Mean of Estimated Parameters:');
disp(mean_theta_hat);
disp('Standard Deviation of Estimated Parameters:');
disp(std_theta_hat);
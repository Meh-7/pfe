% Load the data from the CSV file
data = table2array(readtable('databp.csv'));
% Number of replications
num_replications = 200;

% Store estimated parameters from each replication
theta_hat = zeros(num_replications, 7);
theta0 = zeros(num_replications, 7);

% Loop for multiple replications
for i = 1:num_replications
    % Randomly select starting values for theta0 from a uniform distribution
    theta0(i,:) = rand(1, 7);  % Values between 0 and 1
   
    % Define lower and upper bounds for the parameters
    lb = [0; 0; 0; 0; 0; 0; 0];  
    ub = [2; 1; 1; 2; 1; 1; 2]; 
    
    % Randomly perturb the ranges for generating initial theta0
    lb_perturbed = lb - rand(size(lb)) * 0.1;  % Perturb lower bounds
    ub_perturbed = ub + rand(size(ub)) * 0.1;  % Perturb upper bounds
    
    % Randomly select starting values for theta0 from the perturbed ranges
    theta0(i,:) = lb_perturbed + rand(size(lb_perturbed)) .* (ub_perturbed - lb_perturbed);
    
    % Ensure theta0 remains within bounds
    theta0(i,:) = max(lb.', theta0(i,:));
    theta0(i,:) = min(ub.', theta0(i,:));
    
    % Set fmincon options using optimset
    options = optimset('fmincon');
    options.Display = 'off';
    options.GradObj = 'on'; % Specify that the gradient is provided by the objective function
    
    % Estimate the parameters using fmincon
    [theta_hat(i,:), fval, exitflag, output] = fmincon(@(theta) RbivariateINGARCHnll(theta, data), theta0(i,:), [], [], [], [], lb, ub, [], options);
end
end

% Calculate the mean and standard deviation for each parameter
mean_theta_hat = mean(theta_hat);
std_theta_hat = std(theta_hat);

% Display results

disp('Mean of Estimated Parameters:');
disp(mean_theta_hat);
disp('Standard Deviation of Estimated Parameters:');
disp(std_theta_hat);
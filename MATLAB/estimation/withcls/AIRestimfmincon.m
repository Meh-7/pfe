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
    
    % Generate valid starting values ensuring stationarity
    while true 
        % Generate omega values
        theta0(i,1) = lb(1) + rand() * (ub(1) - lb(1));
        theta0(i,6) = lb(6) + rand() * (ub(6) - lb(6));

        % Generate first element of each row of A and B
        theta0(i,2) = lb(2) + rand() * (ub(2) - lb(2));
        theta0(i,4) = lb(4) + rand() * (ub(4) - lb(4));
        theta0(i,7) = lb(7) + rand() * (ub(7) - lb(7));
        theta0(i,9) = lb(9) + rand() * (ub(9) - lb(9)); 

        % Calculate remaining space for second element in each row to ensure sum < 1
        remaining1A = 1 - theta0(i,2); 
        remaining1B = 1 - theta0(i,4);
        remaining2A = 1 - theta0(i,7); 
        remaining2B = 1 - theta0(i,9);

        % Generate second element of each row using remaining space
        theta0(i,3) = lb(3) + rand() * remaining1A * (ub(3) - lb(3)); 
        theta0(i,5) = lb(5) + rand() * remaining1B * (ub(5) - lb(5)); 
        theta0(i,8) = lb(8) + rand() * remaining2A * (ub(8) - lb(8)); 
        theta0(i,10) = lb(10) + rand() * remaining2B * (ub(10) - lb(10)); 

        % Generate delta
        theta0(i,11) = lb(11) + rand() * (ub(11) - lb(11)); 

        % Check if all row sums are less than 1
        if (theta0(i,2) + theta0(i,3) < 1) && (theta0(i,4) + theta0(i,5) < 1) && ...
           (theta0(i,7) + theta0(i,8) < 1) && (theta0(i,9) + theta0(i,10) < 1)
            break; % Exit loop if starting values are valid
        end
    end
    
    %theta0(i,:) = lb + rand(size(lb)).* (ub - lb);
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
% Load the data from the CSV file
data = table2array(readtable('databp.csv'));

% Define fixed parameters (you need to set these values)
omega1 = 0.3; %mean_theta_hat(1);
alpha11 = 0.2; %mean_theta_hat(2);
beta11 = 0.5; %mean_theta_hat(3);
omega2 = 0.5; %mean_theta_hat(4);
alpha22 = 0.4; %mean_theta_hat(5);
%beta22 = 0.3; %mean_theta_hat(6);
delta = 0.7; %mean_theta_hat(7);

% Create a range of beta22 values
x_values = linspace(0, 1.2, 100);  % Example: 100 values between 0 and 1

% Store the negative log-likelihood values
nll_values = zeros(length(x_values), 1);

% Loop through the beta22 values
for i = 1:length(x_values)
    % Create theta vector with fixed parameters and current beta22
    theta = [omega1; alpha11; beta11; omega2; alpha22; x_values(i); delta];

    % Calculate the negative log-likelihood 
    [nll_values(i), ~, ~] = bivariateINGARCHnll(theta, data); 
end

% Plot the negative log-likelihood against beta22 values
figure;
plot(x_values, nll_values);
xlabel('beta222');
ylabel('Negative Log-Likelihood');
title('Negative Log-Likelihood');
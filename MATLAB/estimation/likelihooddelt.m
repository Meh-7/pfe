% Load the data from the CSV file
data = table2array(readtable('databp2.csv'));

% Define fixed parameters (you need to set these values)
omega1 = 0.3; %mean_theta_hat(1);
alpha11 = 0.2; %mean_theta_hat(2);
beta11 = 0.5; %mean_theta_hat(3);
omega2 = 0.5; %mean_theta_hat(4);
alpha22 = 0.4;%mean_theta_hat(5);
%beta22 = 0.3;%mean_theta_hat(6);
%delta = 0.7; %mean_theta_hat(7);

% Create a range of values
x_values = linspace(0, 1, 100);  % Example: 100 values between -1 and 1
y_values = linspace(0, 1, 100);
% Store the negative log-likelihood values
nll_values = zeros(length(x_values), length(y_values));
% Loop through the values
for i = 1:length(x_values)
    for j = 1:length(y_values)
        % Create theta vector with fixed parameters and current delta
        theta = [omega1; alpha11; beta11 ; omega2; alpha22; x_values(i); y_values(j)];

        % Calculate the negative log-likelihood 
        [nll_values(i,j), ~, ~] = bivariateINGARCHnll(theta, data); 

    end
end

% Plot the negative log-likelihood against delta values
figure;
surf(alpha11_values, omega_values, nll_values);
xlabel('beta22');
ylabel('delta');
zlabel('nll');
title('Negative Log-Likelihood');

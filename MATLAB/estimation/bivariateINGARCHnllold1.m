function [nll, lambda1, lambda2] = bivariateINGARCHnll(theta, data)
    omega1 = theta(1);
    alpha11 = theta(2);
    beta11 = theta(3);
    omega2 = theta(4);
    alpha22 = theta(5);
    beta22 = theta(6);
    delta = theta(7);

    n = size(data, 1);
    lambda1 = zeros(n, 1);
    lambda2 = zeros(n, 1);

     % Initialize lambda values
     lambda1(1) = mean(data(1:3, 1)); % Mean of first 3 data points for series 1 
     lambda2(1) = mean(data(1:3, 2)); % Mean of first 3 data points for series 22);

    % Calculate lambda values for each time step
    for t = 2:n
        lambda1(t) = omega1 + alpha11 * lambda1(t-1) + beta11 * data(t-1, 1);
        lambda2(t) = omega2 + alpha22 * lambda2(t-1) + beta22 * data(t-1, 2);
    end

     % Calculate the log-likelihood function
    c = 1 - exp(-1);
    phi = 1 + delta * (exp(-data(:, 1)) - exp(-c * lambda1)) .* (exp(-data(:, 2)) - exp(-c * lambda2));
    
    
    % Check for NaN or Inf in lambda and phi
    if any(isnan([lambda1; lambda2])) || any(isinf([lambda1; lambda2])) || any(isnan(phi)) || any(isinf(phi))
        nll = 1e10;  % Return a very large value if undefined values are found
        return;
    end

    nll = -sum(data(:, 1) .* log(lambda1) + data(:, 2) .* log(lambda2) - lambda1 - lambda2 + log(phi));
    
    % Add penalty for violating the delta constraint
    for t = 1:n
        constraint_value = (abs(delta) - 1 / ((1 - exp(-c * lambda1(t))) * (1 - exp(-c * lambda2(t)))) ) / abs(delta);
        if constraint_value > 0
            nll = nll + 1e9 * constraint_value;  % Large penalty for violation
        end
    end
end

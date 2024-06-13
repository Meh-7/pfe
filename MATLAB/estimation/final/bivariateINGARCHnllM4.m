function [nll, grad] = bivariateINGARCHnllM4(theta, data)
    omega1 = theta(1);
    alpha11 = theta(2);
    alpha12 = theta(3);
    beta11 = theta(4);
    omega2 = theta(5);
    alpha21 = theta(6);
    alpha22 = theta(7);
    beta22 = theta(8);
    delta = theta(9);

    n = size(data, 1);
    lambda1 = zeros(n, 1);
    lambda2 = zeros(n, 1);
    grad = zeros(size(theta));

    % Initialize lambda values with the first data point
    lambda1(1) = data(1, 1);
    lambda2(1) = data(1, 2);

    % Initialize derivatives of lambda w.r.t. parameters
    d_lambda1_d_omega1 = zeros(n, 1);
    d_lambda1_d_alpha11 = zeros(n, 1);
    d_lambda1_d_alpha12 = zeros(n, 1);
    d_lambda1_d_beta11 = zeros(n, 1);
    d_lambda2_d_omega2 = zeros(n, 1);
    d_lambda2_d_alpha21 = zeros(n, 1);
    d_lambda2_d_alpha22 = zeros(n, 1);
    d_lambda2_d_beta22 = zeros(n, 1);
    
    d_lambda1_d_omega1(1) = 1;
    d_lambda2_d_omega2(1) = 1;


    % Calculate lambda values for each time step
    for t = 2:n
        lambda1(t) = omega1 + alpha11 * lambda1(t-1) + alpha12 * lambda2(t-1) + beta11 * data(t-1, 1);
        lambda2(t) = omega2 + alpha21 * lambda1(t-1) + alpha22 * lambda2(t-1) + beta22 * data(t-1, 2);

        % Calculate derivatives of lambda w.r.t. parameters
        d_lambda1_d_omega1(t) = 1 + alpha11 * d_lambda1_d_omega1(t-1);
        d_lambda1_d_alpha11(t) = lambda1(t-1) + alpha11 * d_lambda1_d_alpha11(t-1);
        d_lambda1_d_alpha12(t) = lambda2(t-1) + alpha11 * d_lambda1_d_alpha12(t-1);
        d_lambda1_d_beta11(t) = data(t-1, 1) + alpha11 * d_lambda1_d_beta11(t-1);
        d_lambda2_d_omega2(t) = 1 + alpha22 * d_lambda2_d_omega2(t-1);
        d_lambda2_d_alpha21(t) = lambda1(t-1) + alpha22 * d_lambda2_d_alpha21(t-1);
        d_lambda2_d_alpha22(t) = lambda2(t-1) + alpha22 * d_lambda2_d_alpha22(t-1);
        d_lambda2_d_beta22(t) = data(t-1, 2) + alpha22 * d_lambda2_d_beta22(t-1);
    end

    % Calculate the log-likelihood function
    c = 1 - exp(-1);
    phi = 1 + delta .* (exp(-data(:, 1)) - exp(-c * lambda1)) .* (exp(-data(:, 2)) - exp(-c * lambda2));

    % Regularization to avoid log(0)
    lambda1 = max(lambda1, eps);
    lambda2 = max(lambda2, eps);
    phi = max(phi, eps);

    % Check for NaN or Inf in lambda and phi
    if any(isnan([lambda1; lambda2])) || any(isinf([lambda1; lambda2])) || any(isnan(phi)) || any(isinf(phi))
        nll = 1e10;  % Return a very large value if undefined values are found
        grad = [];
        return;
    end

    nll = -sum(data(:, 1) .* log(lambda1) + data(:, 2) .* log(lambda2) - lambda1 - lambda2 + log(phi));

    % Calculate derivatives of the log-likelihood function w.r.t. parameters
    for t = 1:n
        % Common terms
        term1 = (data(t, 1) / lambda1(t)) - 1 + ((delta*c*exp(-c*lambda1(t))*(exp(-data(t,2))-exp(-c*lambda2(t))))/phi(t));
        term2 = (data(t, 2) / lambda2(t)) - 1 +((delta*c*exp(-c*lambda2(t))*(exp(-data(t,1))-exp(-c*lambda1(t))))/phi(t));

        % Gradient w.r.t theta_1
        grad(1) = grad(1) + (term1 * d_lambda1_d_omega1(t));
        grad(2) = grad(2) + (term1 * d_lambda1_d_alpha11(t));
        grad(3) = grad(3) + (term1 * d_lambda1_d_alpha12(t));
        grad(4) = grad(4) + (term1 * d_lambda1_d_beta11(t));
        % Gradient w.r.t theta_2
        grad(5) = grad(5) + (term2 * d_lambda2_d_omega2(t));
        grad(6) = grad(6) + (term2 * d_lambda2_d_alpha21(t));
        grad(7) = grad(7) + (term2 * d_lambda2_d_alpha22(t));
        grad(8) = grad(8) + (term2 * d_lambda2_d_beta22(t));

        % Gradient w.r.t delta
        grad(9) = grad(9) + (((exp(-data(t, 1)) - exp(-c * lambda1(t))) * (exp(-data(t, 2)) - exp(-c * lambda2(t)))) / phi(t));
    end

    % Negate the gradient for minimization
    grad = -grad;
end
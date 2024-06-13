function [c,ceq] = phi_constraint(theta, data)
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

    lambda1(1) = omega1 / (1 - alpha11 - beta11);
    lambda2(1) = omega2 / (1 - alpha22 - beta22);

    for t = 2:n
        lambda1(t) = omega1 + alpha11 * lambda1(t-1) + beta11 * data(t-1, 1);
        lambda2(t) = omega2 + alpha22 * lambda2(t-1) + beta22 * data(t-1, 2);
    end

    c = 1 - exp(-1);
    phi = 1 + delta * (exp(-data(:, 1)) - exp(-c * lambda1)) .* (exp(-data(:, 2)) - exp(-c * lambda2));

    % Inequality constraints:  0 < phiL ? |?t | ? phiU
    phiL = 0.1;  % Set your desired lower bound
    phiU = 10;    % Set your desired upper bound
    c = [-phi + phiL; phi - phiU]; % Inequality constraint: c(theta) <= 0 
    ceq = []; 
end

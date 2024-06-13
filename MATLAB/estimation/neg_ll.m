function [ll] = neg_ll(theta, data)
n = size(data,1);
    omega_1 = theta(1);
    alpha_11 = theta(2);
    alpha_12 = theta(3);
    beta_11 = theta(4);
    beta_12 = theta(5);
    
    omega_2 = theta(6);
    alpha_21 = theta(7);
    alpha_22 = theta(8);
    beta_21 = theta(9);
    beta_22 = theta(10);
    delta = theta(11);
    
% Assemble omega, A, B from the destructured values
omega = [omega_1, omega_2];
A = [alpha_11, alpha_12; alpha_21, alpha_22];
B = [beta_11, beta_12; beta_21, beta_22];
% Initialize lambda array 
lambda = zeros(2, n);
% Initialize phi array 
phi = zeros(1, n);
c = 1 - exp(-1);
ll = 0;
 for i = 2:n
        % Update lambda with the recursive formula
        lambda(:, i) = transpose(omega) + A * lambda(:, i - 1) + B * transpose(data(i - 1, :));
        
        % Ensure lambda is non-negative
        lambda(lambda < 0) = 1e-8;  % Avoid negative values
        
        % Calculate phi
        phi(i) = 1 + delta * (exp(-data(i, 1)) - exp(-c * lambda(1, i))) *(exp(-data(i, 2)) - exp(-c * lambda(2, i)));
        %disp(phi(i))
        % Log-likelihood computation
        ll = ll + data(i, 1) * log(lambda(1, i)) + data(i, 2) * log(lambda(2, i)) -lambda(1,i) - lambda(2, i) + log(phi(i));
 end
    
    % Return negative log-likelihood for minimization
    ll = -ll;


end

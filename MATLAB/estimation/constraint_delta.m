function [c, ceq] = constraint_delta(theta, data)
   % No equality constraints
   ceq = []; 
   
   cons = 1 - exp(-1);
   
   % extract lambda
   [~, ~, lambda] = neg_ll_diag(theta, data);
   n = size(lambda, 2);
   
   % Compute the function value for the first column
    lambda1 = lambda(1, 1);
    lambda2 = lambda(2, 1);
    min_value = 1 / ((1 - exp(-cons * lambda1)) * (1 - exp(-cons * lambda2)));
   
   
    % Iterate over the remaining columns of lambda
    for i = 2:n
        % Extract lambda1 and lambda2 from the current column
        lambda1 = lambda(1, i);
        lambda2 = lambda(2, i);

        % Compute the function value
        func_value = 1 / ((1 - exp(-cons * lambda1)) * (1 - exp(-cons * lambda2)));

        % Update the minimum value if the current function value is smaller
        if func_value < min_value
            min_value = func_value;
        end
    end
    %we go out with the smallest value possible for the bounds of delta
    
    % Define bounds for |delta|
    lower_bound = -min_value; 
    upper_bound = min_value;
   
   % Extract delta from theta
    delta = theta(7);
   % Inequality constraints to ensure lower_bound <= |phi| <= upper_bound
    c = [
        lower_bound - abs(delta);  % abs(delta) >= lower_bound
        abs(delta) - upper_bound    % abs(delta) <= upper_bound 
    ]; 
end
% Define the custom output function to display log-likelihood at each step
function stop = display_log_likelihood( state)
    % Get the current iteration number and function value
    iteration = state.iteration;
    neg_ll = state.fval;  % The negative log-likelihood value
    
    % Display the current iteration and the log-likelihood value
    fprintf('Iteration: %d, Negative Log-Likelihood: %.4f\n', iteration, neg_ll);
    
    % No stopping condition (continue optimization)
    stop = false; 
end

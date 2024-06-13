data = table2array(readtable('databp.csv'));
theta0 = [0.225, 0.254, 1, 0.698, 0.349, 0.5, -2+4*rand];
lb = [0; 0; 0; 0; 0; 0; 0];  
ub = [2; 1; 1; 2; 1; 1; 2];
% Set fmincon options using optimset
    options = optimset('fmincon');
    options.Display = 'off';
    options.GradObj = 'on'; % Specify that the gradient is provided by the objective function
    
% Estimate the parameters using fmincon
    [theta_hat, fval, exitflag, output] = fmincon(@(theta) RbivariateINGARCHnll(theta, data), theta0, [], [], [], [], lb, ub, [], options);



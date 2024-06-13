% Define the main folder path
mainFolder = 'M1'; 
n_values = [200; 500; 1000];
N = 200;
% Loop through each subfolder
subfolders = {'size_200', 'size_500', 'size_1000'};
theta0_M1 = [0.6699, 0.1871, 0.2502, 0.0918, 0.1834, 0.4119, 0.1324, 0.2299, 0.3262, 0.2008;
                 0.6018, 0.2164, 0.2488, 0.0942, 0.1574, 0.3421, 0.1031, 0.2865, 0.3471, 0.1935;
                 0.5359, 0.2403, 0.2444, 0.1008, 0.1733, 0.3214, 0.1225, 0.2752, 0.3234, 0.1997];
% Loop over different sample sizes
for i = 1:length(subfolders)
    currentSubfolder = fullfile(mainFolder, subfolders{i});
    n = n_values(i);
    theta_hat = zeros(N-1, size(theta0_M1,2)+1);
    % Loop through each CSV file in the subfolder
    for j = 2:N
         theta0 = [theta0_M1(i,:), -1+2*rand];
        filename = sprintf('databp_%03d.csv', j); 
        filePath = fullfile(currentSubfolder, filename);
        disp(['reading: ', filePath]);
        % Read the CSV file
        % Load the data from the CSV file
        data = table2array(readtable(filePath));
        % Set fmincon options
        options = optimoptions('fmincon', 'Display', 'off', 'GradObj', 'on'); % Suppress output and specify gradient
        lb = [0;0;0;0;0;0;0;0;0;0;-1];  % Lower bounds
        ub = [1;1;1;1;1;1;1;1;1;1;1]; % Upper bounds
        
        % Estimate the parameters using fmincon
        [theta_hat(j-1,:),fval, exitflag, output] = fmincon(@(theta) bivariateINGARCHnllM1(theta, data), theta0, [], [], [], [], lb, ub, [], options);
         % Diagnostic output
        fprintf('Replication %d, n = %d, fval = %.5f, exitflag = %d\n', j-1, n, fval, exitflag);
        fprintf('Estimated Parameters: %s\n', mat2str(theta_hat(j-1,:)));
    end
    % Save theta_hat to a CSV file
    filename = sprintf('theta_hat_M1_%d.csv', n);
    csvwrite(filename, theta_hat);
end

% Define the main folder path
mainFolder = 'M2'; 
n_values = [200; 500; 1000];
N = 200;
% Loop through each subfolder
subfolders = {'size_200', 'size_500', 'size_1000'};

theta0_M2 = [0.363771931, 0.176146651, 0.465447193, 0.663702708, 0.309879223, 0.294166644;
             0.323108229, 0.181653911, 0.485360409, 0.566901453, 0.366247361, 0.290322507;
             0.31728755, 0.191850743, 0.496908549, 0.540703557, 0.377152276, 0.298864612];
% Loop over different sample sizes
for i = 1:length(subfolders)
    currentSubfolder = fullfile(mainFolder, subfolders{i});
    n = n_values(i);
    theta_hat = zeros(N-1, size(theta0_M2,2)+1);
    % Loop through each CSV file in the subfolder
    for j = 2:200
        theta0 = [theta0_M2(i,:), -1+2*rand];
        filename = sprintf('databp_%03d.csv', j); 
        filePath = fullfile(currentSubfolder, filename);
        disp(['reading: ', filePath]);
        % Read the CSV file
        % Load the data from the CSV file
        data = table2array(readtable(filePath));
        % Set fmincon options
        options = optimoptions('fmincon', 'Display', 'off', 'GradObj', 'on'); % Suppress output and specify gradient
        lb = zeros(size(theta0));  % Lower bounds
        ub = [2;1;1;2;1;1;2]; % Upper bounds
        
        % Estimate the parameters using fmincon
        [theta_hat(j-1,:),fval, exitflag, output] = fmincon(@(theta) bivariateINGARCHnllM2(theta, data), theta0, [], [], [], [], lb, ub, [], options);
         % Diagnostic output
        fprintf('Replication %d, n = %d, fval = %.5f, exitflag = %d\n', j-1, n, fval, exitflag);
        fprintf('Estimated Parameters: %s\n', mat2str(theta_hat(j-1,:)));
    end
    % Save theta_hat to a CSV file
    filename = sprintf('theta_hat_M2_%d.csv', n);
    csvwrite(filename, theta_hat);

end

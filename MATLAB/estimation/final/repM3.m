% Define the main folder path
mainFolder = 'M3'; 
n_values = [200; 500; 1000];
N = 200;
% Loop through each subfolder
subfolders = {'size_200', 'size_500', 'size_1000'};

theta0_M3 = [0.333609747, 0.113163145, 0.486662032, 0.451317488, 0.711342914, 0.190078108, -0.04691946, 0.184102729;
             0.297664432, 0.159792818, 0.497712889, 0.372606322, 0.608172202, 0.292373201, 0.02884111, 0.199569528;
             0.314405899, 0.180693143, 0.497698226, 0.351032926, 0.567297391, 0.339279674, 0.046373008, 0.194112692];
% Loop over different sample sizes
for i = 1:length(subfolders)
    currentSubfolder = fullfile(mainFolder, subfolders{i}); 
    n = n_values(i);
    theta_hat = zeros(N-1, size(theta0_M3,2)+1);
    % Loop through each CSV file in the subfolder
    for j = 2:200
        theta0 = [theta0_M3(i,:),(-1+2*rand)];
        filename = sprintf('databp_%03d.csv', j); 
        filePath = fullfile(currentSubfolder, filename);
        disp(['reading: ', filePath]);
        % Read the CSV file
        % Load the data from the CSV file
        data = table2array(readtable(filePath));
        % Set fmincon options
        options = optimoptions('fmincon', 'Display', 'off', 'GradObj', 'on'); % Suppress output and specify gradient
        lb = [0;0;0;0;0;0;0;0;-1];  % Lower bounds
        ub = [1;1;1;1;1;1;1;1;1]; % Upper bounds
        
        % Estimate the parameters using fmincon
        [theta_hat(j-1,:),fval, exitflag, output] = fmincon(@(theta) bivariateINGARCHnllM3(theta, data), theta0, [], [], [], [], lb, ub, [], options);
        % Diagnostic output
        fprintf('Replication %d, n = %d, fval = %.5f, exitflag = %d\n', j-1, n, fval, exitflag);
        fprintf('Estimated Parameters: %s\n', mat2str(theta_hat(j-1,:)));
        
    end
   % Save theta_hat to a CSV file
    filename = sprintf('theta_hat_M3_%d.csv', n);
    csvwrite(filename, theta_hat);
end

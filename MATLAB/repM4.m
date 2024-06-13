% Define the main folder path
mainFolder = 'M4'; 
n_values = [200; 500; 1000];
N = 200;
% Loop through each subfolder
subfolders = {'size_200', 'size_500', 'size_1000'};

theta0_M4 = [0.623313405, 0.139736747, 0.309933421, 0.191621094, 0.644050628, 0.272020279, 0.172661329, 0.289897423;
             0.592485197, 0.162457967, 0.290151599, 0.19964095, 0.457876585, 0.34957847, 0.173476544, 0.296307458;
             0.53910468, 0.18720011, 0.290332972, 0.196256115, 0.361889561, 0.390789927, 0.184106833, 0.299744754];
% Loop over different sample sizes
for i = 1:length(subfolders)
    currentSubfolder = fullfile(mainFolder, subfolders{i});
    n = n_values(i);
    theta_hat = zeros(N-1, size(theta0_M4,2)+1);
    % Loop through each CSV file in the subfolder
    for j = 2:200
        theta0 = [theta0_M4(i,:),-1+2*rand];
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
        [theta_hat(j-1,:),fval, exitflag, output] = fmincon(@(theta) bivariateINGARCHnllM4(theta, data), theta0, [], [], [], [], lb, ub, [], options);
         % Diagnostic output
        fprintf('Replication %d, n = %d, fval = %.5f, exitflag = %d\n', j-1, n, fval, exitflag);
        fprintf('Estimated Parameters: %s\n', mat2str(theta_hat(j-1,:)));
        
    end
    % Save theta_hat to a CSV file
    filename = sprintf('theta_hat_M4_%d.csv', n);
    csvwrite(filename, theta_hat);
end
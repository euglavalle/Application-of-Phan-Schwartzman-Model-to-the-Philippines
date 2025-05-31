delete(gcp('nocreate'));
clear all; clc; close all;
% Open parallel pool (if not already open)
if isempty(gcp('nocreate'))
    parpool('local', 8); % Starts a parallel pool using default settings
end

% Start caffeinate in background and get its PID. Only works on MacOS
[~, pidStr] = system('caffeinate -dims & echo $!');
pid = str2double(strtrim(pidStr)); % Convert PID string to number
disp(['Started caffeinate with PID: ', num2str(pid)]);

% Define parameter ranges
psi_values = 5:0.1:10;  
costy_values = 0.05:0.005:0.1;

% Define target and tolerance
%target_aveby = 35;
%target_aves = 1.6;
%target_aveby = 42;
%target_aves = 2.2;
%tolerance_by = 0.3;
%tolerance_s = 0.1;

% Generate all combinations of parameters
[param_grid_psi, param_grid_costy] = ndgrid(psi_values, costy_values);
num_combinations = numel(param_grid_psi);


% Preallocate results as a struct array
%results(num_combinations) = struct('psi', 0, 'costy', 0, 'aveby', 0, 'aves', 0);
%results = repmat(struct('psi', 0, 'costy', 0, 'aveby', 0, 'aves', 0), num_combinations, 1);
results_matrix = zeros(num_combinations, 4);

% Parallel loop
tic
parfor idx = 1:num_combinations
    disp(['Computing combination ', num2str(idx), ' out of ', num2str(num_combinations)]);
    % Extract parameter values
    psi = param_grid_psi(idx);
    costy = param_grid_costy(idx);
    
    % Run simulation
    out = runme_in_grid_search(costy, psi);
    
    % Store results using indexed assignment
    %results(idx).psi = psi;
    %results(idx).costy = costy;
    %results(idx).aveby = out.base.aveby;
    %results(idx).aves = out.base.aves;
    results_matrix(idx, :) = [psi, costy, out.base.aveby, out.base.aves];
end
toc

% Convert results to table for easier processing
results_table = array2table(results_matrix, 'VariableNames', {'psi', 'costy', 'aveby', 'aves'});

folder_path = 'results_grid_search/PHL';
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

file_path = fullfile(folder_path, 'grid_search_results_phl_scalecur95_18268_shape_5_11328_mugfromEMDAT_newmu_1Cat5storm.csv');
writetable(results_table, file_path);

% Kill caffeinate when done
system(sprintf('kill %d', pid));
disp('Stopped caffeinate.');
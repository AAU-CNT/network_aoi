%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fairness_opf: compare the JFI for the FCFS, HAF, and OPF policies    %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close figures and clear workspace
close all
clearvars

% Define scenario
error = true; % True to evaluate JFI in scenario with errors
plotting = true; % Turn on plots


sources = 'ksources';
errors = 'noerror';
if (error)
    errors = 'error';
end

if (plotting)
    f1 = figure(1);
    hold on
end

fcfs_file = strcat('results_', sources, '_fifo_', errors, '.mat');
opf_file = strcat('results_', sources, '_opf_', errors, '.mat');
haf_file = strcat('results_', sources, '_haf_', errors, '.mat');


% FCFS JFI calculation
load(fcfs_file);
fcfs_jfi = zeros(length(k), length(rho));

for ki = 1 : length(k)
    k = k(ki);
    for j = 1 : length(rho)
        fcfs_jfi(ki, j) = jain_fairness(age_sim_source{ki}(:,j));
    end
end

% OPF JFI calculation
load(opf_file);
opf_jfi = zeros(length(k), length(rho));

for ki = 1 : length(k)
    k = k(ki);
    for j = 1 : length(rho)
        opf_jfi(ki, j) = jain_fairness(age_sim_source{ki}(:,j));
    end
end

% HAF JFI calculation
load(haf_file);
haf_jfi = zeros(length(k), length(rho));

for ki = 1 : length(k)
    k = k(ki);
    for j = 1 : length(rho)
        haf_jfi(ki, j) = jain_fairness(age_sim_source{ki}(:,j));
    end
end


% Save results
save(strcat('jfi_', sources, '_', errors, '.mat'));

if (plotting)
    % Set up color, line and marker scheme
    color= {'b', 'r', 'k', 'c', 'g', 'y','b', 'r', 'k', 'c'};
    markers = {'s', 'o', 'd', 's', 'o', 'd', 's', 'o', 'd', 's'};
    linestyles = {'-', ':', '-.', '-', ':', '-.', '-', ':', '-.', '.'};
    
    for ki = 1 : length(k)
        k = k(ki);
        
        legendfcfs = ['FCFS, K =', num2str(k)];
        legendopf = ['OPF, K =', num2str(k)];
        legendhaf = ['HAF, K =', num2str(k)];
        
        % JFI plot
        figure(1)
        plot(rho, fcfs_jfi(ki, :),  color{ki}, 'linestyle', '-','marker', markers{1}, 'markersize', 8, 'linewidth', 1.5, 'DisplayName', legendfcfs)
        plot(rho, opf_jfi(ki, :),  color{ki}, 'linestyle', '-','marker', markers{2}, 'markersize', 8, 'linewidth', 1.5, 'DisplayName', legendopf);
        plot(rho, haf_jfi(ki, :),  color{ki}, 'linestyle', '-','marker', markers{3}, 'markersize', 8, 'linewidth', 1.5, 'DisplayName', legendhaf);
    end
    
    figure(1)
    legend(); xlabel('\rho'); ylabel('JFI');
    
end


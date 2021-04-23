%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     line_runner: simulate a multisource line network and derive AoI     %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Close figures and clear workspace
close all
clearvars

%% Simulation parameters
rho = 0.05 : 0.05 : 0.95; % Traffic load at the last node
k = 6; % Number of links (and potential sources)
num_packets = 100000; % Number of packets per source
rate = 1; % Rate for all links
last_rate = 0.8; % Rate for the last link
cross = true; % True if traffic comes from all sources
policy = 2; % 0 for FCFS, 1 for OPF, 2 for HAF
plotting = false; % True to turn on plots
saving = false; % True to save the simulation results
epsilon = 0 : 0.05 : 0.1;

%% Initialization
if (plotting)
    f1 = figure(1);
    hold on
    f2 = figure(2);
    hold on
end
% Initialize vectors for general results
age_sim = zeros(length(epsilon), length(rho));
delay_sim = zeros(length(epsilon), length(rho));
age_lower = zeros(length(epsilon), length(rho));
age_upper = zeros(length(epsilon), length(rho));
delay_th = zeros(length(epsilon), length(rho));

% Initialize vectors for source by source analysis
age_sim_source = {length(epsilon)};
delay_sim_source = {length(epsilon)};
paoi_sim_source = {length(epsilon)};
age_lower_source = {length(epsilon)};
age_upper_source = {length(epsilon)};
delay_th_source = {length(epsilon)};
lambda_source = {length(epsilon)};


%% Main simulation cycle
for ki = 1 : length(epsilon)
    mu = [ones(1, k - 1) * rate, last_rate];
    
    age_sim_source{ki} = zeros(k, length(rho));
    delay_sim_source{ki} = zeros(k, length(rho));
    age_lower_source{ki} = zeros(k, length(rho));
    age_upper_source{ki} = zeros(k, length(rho));
    delay_th_source{ki} = zeros(k, length(rho));
    lambda_source{ki} = zeros(k, length(rho));
    
    for j = 1 : length(rho)        
        % Set up arrival rates
        if (cross)
            lambda = ones(1, k) * rho(j) * mu(end) / k;
        else
            lambda = [rho(j) * mu(end), zeros(1, k-1)];
        end
        error = ones(1, k) * epsilon(ki);
        error
        rho(j)

        % Simulation main function
        [av_age_sim, av_age_lower, av_age_upper, system_delay_sim, system_delay_theoretical] = simulate_line(k, lambda, k * ones(1, k + 1), mu, error, num_packets, policy);
        
        % Results registration (by source)
        age_sim_source{ki}(:, j) = av_age_sim;
        delay_sim_source{ki}(:, j) = system_delay_sim;
        age_lower_source{ki}(:, j) = av_age_lower;
        age_upper_source{ki}(:, j) = av_age_upper;
        delay_th_source{ki}(:, j) = system_delay_theoretical;
        lambda_source{ki}(:, j) = lambda;
        
        
        % Compute averages and register results
        age_sim(ki, j) = sum(av_age_sim .* lambda) / sum(lambda);
        age_lower(ki, j) = sum(av_age_lower .* lambda) / sum(lambda);
        age_upper(ki, j) = sum(av_age_upper .* lambda) / sum(lambda);
        delay_sim(ki, j) = sum(system_delay_sim .* lambda) / sum(lambda);
        delay_th(ki, j) = sum(system_delay_theoretical .* lambda) / sum(lambda);
    end
    
end

%% Plot AoI and system delay
if (plotting)
    % Set up color, line and marker scheme
    color= {'b', 'r', 'k', 'c', 'g', 'y','b', 'r', 'k', 'c'};
    markers = {'s', 'o', 'd', 's', 'o', 'd', 's', 'o', 'd', 's'};
    linestyles = {'-', ':', '-.', '-', ':', '-.', '-', ':', '-.', '.'};
    
    for ki = 1 : length(epsilon)
        e = epsilon(ki);
        
        legendsim = ['Simulation, \epsilon =', num2str(e)];
        legendthe = ['Bounds, \epsilon =', num2str(e)];
        
        % Average AoI plot
        figure(1)
        plot(rho, age_sim(ki, :), color{ki}, 'linestyle', 'none','marker', markers{ki}, 'markersize', 8, 'linewidth', 1.5, 'DisplayName', legendsim)
        plot(rho, age_lower(ki, :),  color{ki}, 'linestyle', '-', 'linewidth', 1.5, 'DisplayName', legendthe)
        h = plot(rho, age_upper(ki, :),  color{ki}, 'linestyle', '-', 'linewidth', 1.5, 'DisplayName', legendthe);
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
        
        % System delay plot
        figure(2)
        plot(rho, delay_sim(ki, :), color{ki}, 'linestyle', 'none','marker', markers{ki}, 'markersize', 8, 'linewidth', 1.5, 'DisplayName', legendsim)
        legendthe = ['\epsilon =', num2str(e)];
        plot(rho, delay_th(ki, :),  color{ki}, 'linestyle', '-', 'linewidth', 1.5, 'DisplayName', legendthe)
        
    end
  
    figure(1)
    legend();
    xlabel('\rho');
    ylabel('Average AoI');
    
    figure(2)
    legend();
    xlabel('\rho');
    ylabel('Average system delay');
end

%% Save results
if (saving)
    sources = 'singlesource';
    if (cross)
        sources = 'ksources';
    end
    policy_text = '_fifo_';
    error = 'noerror';
    if (sum(epsilon) > 0)
        error = 'error';
    end
    if (policy == 1)
        policy_text = '_opf_';
    end
    if (policy == 2)
        policy_text = '_haf_';
    end
    save(strcat('results_', sources, policy_text, error, '.mat'));
end

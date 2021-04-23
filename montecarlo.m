function [av_age_sim, delay_sim] = montecarlo(k, lambda, mu_isl, mu_downlink, epsilon_isl, epsilon_uplink, epsilon_downlink, num_packets, policy, rho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    function simulate: main simulation loop, computes theoretical and    %
%          Monte Carlo delay and AoI in the Poisson tandem queue          %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% k : number of satellites                                                %
% lambda (1*k): ground source arrival rate vector                         %
% mu_isl (1*(k-1)): inter-satellite link service rate vector              %
% mu_downlink: satellite to ground link service rate                      %
% epsilon_isl (1*(k-1)): inter-satellite link service rate vector         %
% epsilon_uplink (1*k): ground to satellite link service rate vector      %
% epsilon_downlink: satellite to ground link service rate                 %
% num_packets: number of packets to consider in Monte Carlo simulation    %
% policy: 0 for FCFS, 1 for OPF, 2 for HAF                                %
%                                                                         %
% av_age_sim (1*k): Monte Carlo average AoI for each source               %
% delay_sim (1*k): Monte Carlo average system delay for each source       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

origin_times = ones(k, num_packets) * 1E100;

if (rho > 0)
    load('departure_times_aloha2.mat', 'pcktDepartureTime_Aloha_sim');
    rho = round(rho);
    for i = 1 : k
        if (lambda(i) > 0)
            origin = (pcktDepartureTime_Aloha_sim(rho, num_packets * (i - 1) + 1 : num_packets * i) * k) - 1;
            if (i > 1)
                origin = origin - max(pcktDepartureTime_Aloha_sim(rho, 1 : num_packets * (i - 1)) * k) ;
            end
            origin(origin<0) = 1E100;
            origin = sort(origin);
            interarr = diff(origin(origin<1E100));
            mult = 1 / lambda(i) / mean(interarr);
            origin(origin<1E100) = origin(origin<1E100) * mult;
            interarr = diff(origin(origin<1E100));
            origin_times(i, :) = origin;
        end
    end
else
    % Compute origin times at the sources
    for j = 1 : k
        if (lambda(j) > 0)
            origin_times(j, :) = cumsum((exprnd(1 / lambda(j), 1, num_packets)));
        end
    end
end


% Initialization of packet matrices
arrival_times = -ones(k, num_packets, k);

min_index = 1000; % Parameter to avoid transitory behavior

% Results
av_age_sim = zeros(1, k);
delay_sim = zeros(1, k);



exit = (ones(1, k) > 0);

% Compute arrival time at each satellite
for j = 1 : k
    % Arrivals from previous satellites
    if (j > 1)
        arrival_times(1 : j - 1, :, j) = multi_source_queue(arrival_times(1 : j - 1, :, j - 1), origin_times(1 : j - 1, :), exit, mu_isl(j - 1), epsilon_isl(j - 1), policy);
    end
    % Arrivals from ground
    arrival_times(j, :, j) = origin_times(j, :);
end
% Compute final departure times from the system (i.e., arrival at sink)
departure_times = multi_source_queue(arrival_times(:, :, k), origin_times, exit, mu_downlink, epsilon_downlink, policy);
departure_times(departure_times == 1E100) = -1;
maxT = min(max(departure_times'));

% Compute average AoI and system delay
for j = 1 : k
    if (lambda(j) > 0)
        [av_age_sim(1, j), delay_sim(1, j)] = average_age(origin_times(j, :), departure_times(j, :), maxT, min_index);
    end
end

end


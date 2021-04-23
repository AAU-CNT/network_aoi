function [av_age_sim, av_age_lower, av_age_upper, delay_sim, delay_theory] = simulate_line(k, lambda, exit, mu, epsilon, num_packets, policy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function simulate_line: main simulation loop, compute theoretical and  %
%          Monte Carlo delay and AoI in the Poisson line network          %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% k: number of links in the network                                       %
% lambda (1*k): source arrival rate vector                                %
% exit (1*k): exit node for each source                                   %
% mu (1*k):  link service rate vector                                     %
% epsilon (1*k): link error rate vector                                   %
% num_packets: number of packets to consider in Monte Carlo simulation    %
% policy: 0 for FCFS, 1 for OPF, 2 for MAF                                %
%                                                                         %
% av_age_sim (1*k): Monte Carlo average AoI for each source               %
% av_age_lower (1*k): theoretical average AoI lower bound for each source %
% av_age_upper (1*k): theoretical average AoI upper bound for each source %
% delay_sim (1*k): Monte Carlo average system delay for each source       %
% delay_theory (1*k): theoretical average system delay for each source    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization of packet matrices
origin_times = ones(k, num_packets) * 1E100;
arrival_times = -ones(k, num_packets, k);

min_index = 1000; % Parameter to avoid transitory behavior

% Results
av_age_sim = zeros(1, k);
av_age_lower = zeros(1, k);
av_age_upper = zeros(1, k);
delay_sim = zeros(1, k);
delay_theory = zeros(1, k);


% Compute origin times at the sources
for j = 1 : k
    if (lambda(j) > 0)
        origin_times(j, :) = cumsum((exprnd(1 / lambda(j), 1, num_packets)));
    end
end


% Compute arrival time at each node
for j = 1 : k
    % Arrivals from previous links
    if (j > 1)
        arrival_times(1 : j - 1, :, j) = multi_source_queue(arrival_times(1 : j - 1, :, j - 1), origin_times(1 : j - 1, :), (exit(1 : j - 1) >= j - 1), mu(j - 1), epsilon(j - 1), policy);
    end
    % Arrivals from new source
    arrival_times(j, :, j) = origin_times(j, :);
end
% Compute final departure times from the system (i.e., arrival at sink)
departure_times = multi_source_queue(arrival_times(:, :, k), origin_times, (exit >= k), mu(k), epsilon(k), policy);
maxtimes = max(origin_times');
maxT = min(maxtimes(maxtimes > min_index));

% Compute average AoI and system delay
for j = 1 : k
    if (lambda(j) > 0 && exit(j) == k)
        [av_age_sim(1, j), delay_sim(1, j)] = average_age(origin_times(j, :), departure_times(j, :), maxT, min_index);
    end
    if (lambda(j) > 0 && exit(j) < k)
        [av_age_sim(1, j), delay_sim(1, j)] = average_age(origin_times(j, :), arrival_times(j, :, exit(j)), maxT, min_index);
    end
end

% Compute lambda vector considering errors
lambdae = zeros(1, k);
lambdae(1) = lambda(1);
for i = 2 : k
    % Consider sources upstream
    for j = 1 : i - 1
        if (i <= exit(j))
            lambdae(i) = lambdae(i) + lambda(j) * prod(1 - epsilon(j : i - 1));
        end
    end
    lambdae(i) = lambdae(i) + lambda(i);
end

% Compute alpha for each node
alpha = zeros(1, k);
for i = 1 : k
    alpha(i) = mu(i) - lambdae(i);
end

% Compute theoretical AoI for each node
for j = 1 : k
    if (lambda(j) > 0)
        % Mu and alpha vectors for source j
        muj = mu(j : min(exit(j), k));
        alphaj = alpha(j : exit(j));
        % Compute success probability for source j
        p_succ = prod(1 - epsilon(j : k));
        % System delay: Hypoexponential average
        delay_theory(j) = sum(1 ./ alphaj);
        av_age_upper(j) = (1 / lambda(j) + delay_theory(j)) / p_succ;
        if (policy > 0)
            % The lower bound is computed using E[WY] for OPF and MAF
            av_age_lower(j) = lambda(j) * opf_WY(alphaj, [lambda(j), lambdae(j : k)], muj);
        else
            % The lower bound is computed using E[WY] for the FCFS policy
            av_age_lower(j) = lambda(j) * hypoexp_WY(alphaj, lambda(j), muj(1 : end - 1));
        end
        av_age_lower(j) = av_age_lower(j) + 1 / (lambda(j) * p_succ);
        av_age_lower(j) = av_age_lower(j) + sum(1 ./ muj);
        av_age_lower(j) = av_age_lower(j) + (1 - p_succ) * (1 / muj(length(muj)) - 2 / lambda(j));
    end
end

end


function [av_age_sim, av_age_lower, av_age_upper, delay_sim, delay_theory, tail_sim, tail_upper] = simulate_dumbbell(k, lambda, common, mu, epsilon, num_packets, policy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    function simulate: main simulation loop, computes theoretical and    %
%         Monte Carlo delay and AoI in the Poisson dumbbell queue         %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% k : number of satellites                                                %
% lambda (1*n): source arrival rate vector                                %
% common: bottleneck link                                                 %
% mu (1*k): link service rate vector                                      %
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

% Number of sources
N = length(lambda);

% Initialization of packet matrices
origin_times = ones(N, num_packets) * 1E100;
arrival_times = -ones(N, num_packets, k + 1);

min_index = 1000; % Parameter to avoid transitory behavior

% Results
av_age_sim = zeros(1, N);
av_age_lower = zeros(1, N);
av_age_upper = zeros(1, N);
delay_sim = zeros(1, N);
delay_theory = zeros(1, N);

% Compute origin times at the sources
for n = 1 : N
    origin_times(n, :) = cumsum((exprnd(1 / lambda(n), 1, num_packets)));
    arrival_times(n, :, 1) = origin_times(n, :);
end


% Compute arrival time at each link before the bottleneck
for n = 1 : N
    for j = 2 : common
        % Arrivals from previous links
        arrival_times(n, :, j) = multi_source_queue(arrival_times(n, :, j - 1), origin_times(n, :), 1, mu(j - 1), epsilon(j - 1), policy);
    end
end
% Compute arrival time just after the bottleneck
arrival_times(:, :, common + 1) = multi_source_queue(arrival_times(:, :, common), origin_times, ones(1, N), mu(common), epsilon(common), policy);
% Compute arrival time at each link after the bottleneck
for n = 1 : N
    for j = common + 2 : k + 1
        arrival_times(n, :, j) = multi_source_queue(arrival_times(n, :, j - 1), origin_times(n, :), 1, mu(j - 1), epsilon(j - 1), policy);
    end
end

% Compute final departure times from the system (i.e., arrival at sink)
departure_times = squeeze(arrival_times(:, :, k + 1));
maxtimes = max(origin_times');
maxT = min(maxtimes(maxtimes > min_index));

% Compute average AoI and system delay
for j = 1 : N
    [av_age_sim(1, j), delay_sim(1, j)] = average_age(origin_times(j, :), departure_times(j, :), maxT, min_index);
end

% Compute lambda vector considering errors
lambdae = zeros(n, k);
lambdae(:, 1) = lambda;
for i = 2 : k
    lambdae(:, i) = lambdae(:, i - 1) * (1 - epsilon(i));
end
lambdae(:, common) = sum(lambdae(:, common));

% Compute alpha for each node
alpha = zeros(n, k);
for i = 1 : k
    alpha(:, i) = mu(i) - lambdae(:, i);
end

% Compute theoretical AoI for each source
for n = 1 : N
    if (lambda(n) > 0)
        % Mu and alpha vectors for source n
        alphan = alpha(n, :);
        % Compute success probability for source n
        p_succ = prod(1 - epsilon);
        % System delay: Hypoexponential average
        delay_theory(n) = sum(1 ./ alphan);
        av_age_upper(n) = (1 / lambda(n) + delay_theory(n)) / p_succ;
        if (policy > 0)
            % The lower bound is computed using E[WY] for the OPF and MAF
            % policied
            av_age_lower(n) = lambda(n) * opf_WY(alphan, lambdae(1 : k), mu);
        else
            % The lower bound is computed using E[WY] for the FCFS policy
            av_age_lower(n) = lambda(n) * hypoexp_WY(alphan, lambda(n), mu(1 : end - 1));
        end
        av_age_lower(n) = av_age_lower(n) + 1 / (lambda(n) * p_succ);
        av_age_lower(n) = av_age_lower(n) + sum(1 ./ mu(1 : end - 1));
        av_age_lower(n) = av_age_lower(n) + (1 - p_succ) * (1 / mu(k) - 2 / lambda(n));
    end
end

end


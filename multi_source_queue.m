function [departure_times] = multi_source_queue(arrival_times, origin_times, exit, mu, epsilon, policy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function multi_source_queue: computes the departure times for a queue  %
%    with link errors and multiple sources, using FIFO or OPF policies    %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% arrival_times (k*N): arrival time for each packet from each source (-1  %
%                      indicates a dropped packet)                        %
% origin_times (k*N): generation time for each packet from each source    %
% exit (1*k): exit bool for each source                                   %
% mu: departure rate for the queue                                        %
% epsilon: error rate of the link                                         %
% policy: 0 for FCFS, 1 for OPF, 2 for HAF                                %
%                                                                         %
% departure_times (k*N): departure time for each packet from each source  %
%                        (-1 indicates a dropped packet)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization of support variables and output vector
time = 0;
k = size(arrival_times, 1);
num_packets = size(arrival_times, 2);
departure_times = ones(k, num_packets) * 1E100;
% Start from the first packet from each source
next_index = ones(1, k);
% Initialization of priority vector

for j = 1 : k
    if (~exit(j))
        arrival_times(j, :) = -1;
    end
end

if (policy == 0)
    % FCFS: priority is the same as the arrival time
    priority = arrival_times(:, 1);
end

if (policy > 0)
    % In OPF/MAF, priority is different, but the first packet to arrive is
    % the first to be served
    [time, first] = min(arrival_times(:, 1));
    priority = ones(1, k) * 1E100;
    priority(first) = time;
end

while (min(next_index) <= num_packets)
    % Find the first packet in the queue
    [~, next_source] = min(priority);
    % If multiple packets have the same priority, choose one at random
    if (length(next_source) > 1)
        next_source = next_source(randi(length(next_source)));
    end
    % End the cycle: all packets have been delivered
    if(next_index(next_source) > num_packets)
        break;
    end
    % Deliver the packet and update time step
    time = max(time, arrival_times(next_source, next_index(next_source))) + exprnd(1 / mu);
    % Insert link error
    if (rand > epsilon)
        departure_times(next_source, next_index(next_source)) = time;
    else
        departure_times(next_source, next_index(next_source)) = -1;
    end
    % Increase index for the source of the delivered packet
    next_index(next_source) = next_index(next_source) + 1;
    % Discard all packets with errors
    while (next_index(next_source) <= num_packets && arrival_times(next_source, next_index(next_source)) < 0)
        departure_times(next_source, next_index(next_source)) = -1;
        next_index(next_source) = next_index(next_source) + 1;
    end
    % Update priority vector
    if (policy > 0)
        % Keep track of next arrival
        next_arrival = 1E100;
        arr_index = 0;
        % Iterate over sources
        for source = 1 : k
            % Check if all packets from the source have been delivered
            if (next_index(source) > num_packets)
                priority(source) = 1E100;
            else
                % If empty is true, there are no queued packets from
                % this source
                empty = time < arrival_times(source, next_index(source));
                if (empty)
                    priority(source) = 1E100;
                    % Update next arrival
                    if (arrival_times(source, next_index(source)) < next_arrival)
                        next_arrival = arrival_times(source, next_index(source));
                        arr_index = source;
                    end
                else
                    if (policy == 1)
                        % In OPF, priority is the origin time
                        priority(source) = origin_times(source, next_index(source));
                    else
                        % In MAF, priority is the AoI (ie, origin time of
                        % the previous packet)
                        if (next_index(source) > 1)
                            priority(source) = origin_times(source, next_index(source) - 1);
                        else
                            priority(source) = 0;
                        end
                    end
                end
            end
        end
        % The system is empty, the first served packet is the next one
        % to arrive
        if (min(priority) == 1E100 && arr_index > 0)
            priority(arr_index) = next_arrival;
        end
    else
        % FCFS: priority is the same as arrival time
        priority(next_source) = 1E100;
        if (next_index(next_source) <= num_packets)
            priority(next_source) = arrival_times(next_source, next_index(next_source));
        end
    end
end


end


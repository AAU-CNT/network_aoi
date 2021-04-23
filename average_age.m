function [av_age, av_delay] = average_age(origin_times, departure_times, maxT, min_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function average_age: computes the average AoI and system delay     %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% origin_times (1*N): generation time for each packet                     %
% departure_times (1*N): system departure time for each packet (a -1      %
%                        value indicates the packet was lost)             %
% maxT: maximum time to evaluate                                          %
% min_index: minimum index to evaluate (to remove initial transition)     %
%                                                                         %
% av_age: average AoI                                                     %
% av_delay: average system delay                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
av_age = 0;
av_delay = 0;
t = 0;
n = 0;

% Main evaluation loop
for i = min_index : length(origin_times) - 1
    % Check for errors in packet reception
    if (departure_times(i) > 0)
        % Find next correctly received packet
        j = 1;
        while (i + j < length(origin_times) && departure_times(i + j) < 0)
            j = j + 1;
        end
        % Check maximum time condition
        if (departure_times(i + j) < maxT)
            % Compute average age in the time between two packets
            age = (departure_times(i + j) + departure_times(i) - 2 * origin_times(i)) / 2;
            interdep = departure_times(i + j) - departure_times(i);
            av_age = av_age + age * interdep;
            % Add delay for the packet to average
            av_delay = av_delay + departure_times(i) - origin_times(i);
            % Increase total time and number of samples
            t = t + interdep;
            n = n + 1;
        else
            break
        end
    end
end
% Division by total time/number of samples to get average
av_delay = av_delay / n;
av_age = av_age / t;

end
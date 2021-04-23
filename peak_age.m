function [age] = peak_age(origin_times, departure_times, maxT, min_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 function peak_age: compute the peak AoI                 %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% origin_times (1*N): generation time for each packet                     %
% departure_times (1*N): system departure time for each packet (a -1      %
%                        value indicates the packet was lost)             %
% maxT: maximum time to evaluate                                          %
% min_index: minimum index to evaluate (to remove initial transition)     %
%                                                                         %
% age: PAoI vector                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
age = [];

% Main evaluation loop
for i = min_index : length(origin_times) - 1
    % Check for errors in packet reception
    if (departure_times(i) > 0)
        % Find next correctly received packet
        j = 1;
        while (i + j <= length(origin_times) && departure_times(i + j) < 0)
            j = j + 1;
        end
        % Check maximum time condition
        if (departure_times(i + j) < maxT)
            % Compute average age in the time between two packets
            next_age = (departure_times(i + j) + departure_times(i) - 2 * origin_times(i)) / 2;
            age = [age, next_age];
        else
            break
        end
    end
end

end
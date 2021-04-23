function [sub] = hypoexp_coeff_sub(tot, alpha, k, i, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function hypoexp_coeff_sub: get Hypoexponential sub-coefficient     %
% (recursive function iterating over possible vectors that sum to k(i)-j) %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% tot: the total sum left for the remaining elements of m                 %
% alpha (1*L): the system time rate parameter vector (unique)             %
% k (1*L): the occurrence vector for alpha                                %
% i: the selected alpha for the coefficient                               %
% m (1*(M-L)): the parameter vector for the sub-coefficient               %
%                                                                         %
% sub: the sub-coefficient for the Hypoexponential                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the sub-coefficient
sub = 0;

% Recursive function
if (length(m) < length(alpha) - 1)
    % Check if the next index is i: m(i) is always 0
    if (length(m) == i - 1)
        sub = sub + hypoexp_coeff_sub(tot, alpha, k, i, [m 0]);
    else
        % Iterate over all possible values of the next m
        for next_m = 0 : tot
            % Get recursive sub-coefficient
            sub = sub + hypoexp_coeff_sub(tot - next_m, alpha, k, i, [m next_m]);
        end
    end
% Base case for the recursion
else
    % Check if i is the last index
    if (tot > 0 && length(m) == i - 1)
        % The vector does not sum to k(i)-j: discard case
        sub = 0;
    else
        % Complete m so it sums to k(i)-j
        m = [m tot];
        sub = 1;
        % Compute sub-coefficient for the given m vector
        for j = 1 : length(alpha)
            if (j ~= i)
                sub = sub * nchoosek(k(j) + m(j) - 1, m(j)) / ((alpha(j) - alpha(i)) ^ (k(j) + m(j)));
            end
        end
    end
end

end


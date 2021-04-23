function [P] = hypoexp_cdf(x, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hypoexp_cdf: compute the CDF of a Hypoexponential distribution %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% x (1*n): vector of points                                               %
% lambda (1*k): parameter vector of the distribution                      %
%                                                                         %
% p (1*n): PDF of the hypoexponential in the points of x                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


P = zeros(size(x));
[lambda, n] = compact_vec(lambda);


% Hypoexponential sum
for i = 1 : length(lambda)
    for j = 1 : n(i)
        gamma = hypoexp_coeff(lambda, n, i, j);
        % Iterate over points
        for point = 1 : length(x)
            P(point) = P(point) + gamma * 1 / lambda(i) ^ j;
            for l = 0 : j - 1
                P(point) = P(point) - gamma * exp(-lambda(i) * x(point)) * x(point) ^ l / factorial(l) / lambda(i) ^ (j - l);
            end
        end
    end
end

end


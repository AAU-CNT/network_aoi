function [p] = hypoexp_pdf(x, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hypoexp_pdf: compute the PDF of a Hypoexponential distribution %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% x (1*n): vector of points                                               %
% lambda (1*k): parameter vector of the distribution                      %
%                                                                         %
% p (1*n): PDF of the hypoexponential in the points of x                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p = zeros(size(x));
[lambda, n] = compact_vec(lambda);


% Hypoexponential sum
for i = 1 : length(lambda)
    for j = 1 : n(i)
        gamma = hypoexp_coeff(lambda, n, i, j);
        % Iterate over points
        for point = 1 : length(x)
            p(point) = p(point) + gamma * exp(-lambda(i) * x(point)) * x(point) ^ (j - 1) / factorial(j - 1);
        end
    end
end

end


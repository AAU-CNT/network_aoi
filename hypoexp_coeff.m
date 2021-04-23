function [coeff] = hypoexp_coeff(alpha, k, i, j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function hypoexp_coeff: get Hypoexponential coefficient         %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% alpha (1*M): the system time rate parameter vector (unique)             %
% k (1*M): the occurrence vector for alpha                                %
% i: the selected alpha for the coefficient                               %
% j: the exponent for the coefficient (between 1 and k(i))                %
%                                                                         %
% compact (1*M): the vector with only unique elements                     %
% occurrences (1*M): the vector with occurrence numbers for each element  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize coefficient
coeff = (-1) ^ (k(i) - j);

% Multiply by product of all alpha values
for l = 1 : length(alpha)
    coeff = coeff * (alpha(l) ^ k(l));
end

% Get sub-coefficient
coeff = coeff * hypoexp_coeff_sub(k(i) - j, alpha, k, i, []);

end


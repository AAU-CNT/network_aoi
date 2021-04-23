function [wy] = hypoexp_WY(alpha, lambda, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           function hypoexp_WY: compute E[WY] in the FCFS case           %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% alpha (1*N): system time rate parameter vector                          %
% lambda: arrival rate parameter                                          %
% mu (1*N): service rate parameter vector                                 %
%                                                                         %
% wy: E[WY], used to compute average AoI                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wy=0;

% Determine compact vectors
[alpha_c, ka] = compact_vec(alpha);
[mu_c, km] = compact_vec(mu);
Katot = length(alpha_c);
Kmtot = length(mu_c);

% Main computation loop
for i = 1 : Katot
    for j = 1 : ka(i)
        for l = 0 : j - 1
            for m = 0 : l
                for n = 1 : Kmtot
                    for o = 1 : km(n)
                        % Compute Hypoexponential coefficients
                        gamma = hypoexp_coeff(alpha_c, ka, i, j);
                        delta = hypoexp_coeff(mu_c, km, n, o);
                        % Compute contribution to overall AoI
                        num =  gamma * delta * lambda * (j - l) * (l - m + 1) * factorial(m + o - 1);
                        den = (alpha_c(i) ^ (j - l + 1)) * ((alpha_c(i) + lambda) ^ (l - m + 2)) * ((alpha_c(i) + mu_c(n)) ^ (m + o)) * factorial(m) * factorial(o - 1);
                        wy = wy + num / den;
                    end
                end
            end
        end
    end
end
end


function [wy] = opf_WY(alpha, lambda, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           function opf_WY: computes E[WY] in the OPF/MAF case           %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% alpha (1*N): system time rate parameter vector                          %
% lambda (1*N): arrival rate parameter vector                             %
% mu (1*N): service rate parameter vector                                 %
%                                                                         %
% wy: E[WY], used to compute average AoI                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute AoI contribution for first node
wy = lambda(1) / (alpha(1) * (mu(1) ^ 2));

% Compute AoI contribution for later nodes
for j = 2 : length(mu)
    rho = lambda(j) / mu(j);
    if(abs(alpha(j) - mu(j - 1)) < 1E-6)
        wy = wy + 1;
    else
        wy = wy + (1 - rho) * rho * mu(j - 1) / (lambda(j) * (alpha(j) + mu(j - 1)));
    end
end


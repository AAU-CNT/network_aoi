function [jfi] = jain_fairness(age)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function jain_fairness: compute the average AoI Jain Fairness Index   %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% age (1*K): average age vector                                           %
%                                                                         %
% jfi: the Jain Fairness Index between sources                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jfi = (sum(age) ^ 2) / (length(age) * sum(age .^ 2));

end


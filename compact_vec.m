function [compact, occurrences] = compact_vec(vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      function compact_vec: generates a compact version of a vector      %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% vector (1*N): the original vector                                       %
%                                                                         %
% compact (1*M): the vector with only unique elements                     %
% occurrences (1*M): the vector with occurrence numbers for each element  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find unique elements
compact = unique(vector);
if (length(compact) == length(vector))
    occurrences= ones(1, length(vector));
else
    % Count occurrences
    occurrences=zeros(1, length(compact));
    for i = 1 : length(compact)
        occurrences(i) = sum(vector == compact(i));
    end
end

end


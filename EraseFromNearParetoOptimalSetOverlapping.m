function [NearParetoOptimalSet] = EraseFromNearParetoOptimalSetOverlapping(NearParetoOptimalSet)
% EraseFromNearParetoOptimalSetOverlapping
% Removes Pareto entries marked as Dominated = 'T'.

    if isempty(NearParetoOptimalSet)
        return;
    end

    keepIdx = true(1, length(NearParetoOptimalSet));

    for i = 1:length(NearParetoOptimalSet)

        if isfield(NearParetoOptimalSet(i), 'Dominated') && ...
                strcmp(NearParetoOptimalSet(i).Dominated, 'T')

            keepIdx(i) = false;
        end
    end

    NearParetoOptimalSet = NearParetoOptimalSet(keepIdx);

    % Reset Dominated flags for safety
    for i = 1:length(NearParetoOptimalSet)
        NearParetoOptimalSet(i).Dominated = '';
    end
end
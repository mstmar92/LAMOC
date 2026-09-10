function [NearParetoOptimalSet] = UpdateNearParetoOptimalSetOverlapping( ...
    NearParetoOptimalSet, Population, PopSize, ObjectiveDimension)
% UpdateNearParetoOptimalSetOverlapping
% Updates the near Pareto-optimal set for minimization objectives.
%
% Obj(1): Conductance loss       minimize
% Obj(2): GO coherence loss      minimize

    if nargin < 3 || isempty(PopSize)
        PopSize = numel(Population);
    end

    if PopSize > numel(Population)
        error('PopSize is larger than the actual population size.');
    end

    for i = 1:PopSize

        if ~isfield(Population(i), 'Obj') || numel(Population(i).Obj) < ObjectiveDimension
            error('Population(%d) does not contain valid Obj.', i);
        end

        if any(isnan(Population(i).Obj)) || any(isinf(Population(i).Obj))
            continue;
        end

        Flag = 1;

        for nd = 1:length(NearParetoOptimalSet)

            % Avoid duplicate community structures
            if isequal(NearParetoOptimalSet(nd).CmplxID, Population(i).CommsNodes)
                NearParetoOptimalSet(nd).Dominated = '';
                Flag = 0;
                break;
            end

            % If an existing Pareto solution dominates this population member,
            % do not add the current individual.
            if ObjectiveDimension == 2

                oldObj = NearParetoOptimalSet(nd).Obj;
                newObj = Population(i).Obj;

                oldDominatesNew = ...
                    oldObj(1) <= newObj(1) && ...
                    oldObj(2) <= newObj(2) && ...
                    (oldObj(1) < newObj(1) || oldObj(2) < newObj(2));

                if oldDominatesNew
                    NearParetoOptimalSet(nd).Dominated = '';
                    Flag = 0;
                    break;
                end
            else
                error('Only ObjectiveDimension = 2 is currently supported.');
            end

            % Mark existing solutions dominated by the new solution
            NearParetoOptimalSet = Dominate( ...
                Population(i), NearParetoOptimalSet, nd, ObjectiveDimension);
        end

        if Flag

            NearParetoOptimalSet = EraseFromNearParetoOptimalSetOverlapping( ...
                NearParetoOptimalSet);

            L = length(NearParetoOptimalSet);

            NearParetoOptimalSet(L+1).Chromosome = Population(i).Chromosome;
            NearParetoOptimalSet(L+1).Dominated  = '';
            NearParetoOptimalSet(L+1).CmplxID    = Population(i).CommsNodes;
            NearParetoOptimalSet(L+1).Obj        = Population(i).Obj;

            NearParetoOptimalSet(L+1).Obj1_Conductance = Population(i).Obj(1);
            NearParetoOptimalSet(L+1).Obj2_GOCoherenceLoss = Population(i).Obj(2);

            if isfield(Population(i), 'Kv')
                NearParetoOptimalSet(L+1).Kv = Population(i).Kv;
            else
                NearParetoOptimalSet(L+1).Kv = NaN;
            end

            if isfield(Population(i), 'avgGO')
                NearParetoOptimalSet(L+1).avgGO = Population(i).avgGO;
            else
                NearParetoOptimalSet(L+1).avgGO = NaN;
            end
        end
    end
end
function [IdealPoint, IndivPoint] = UpdateReferenceOverlapping(IdealPoint, IndivPoint, ...
    Population, PopulationSize)
% UpdateReferenceOverlapping
% Updates the ideal point for MOEA/D under minimization.
%
% Obj(1): Conductance          minimize
% Obj(2): GO coherence loss    minimize

    if nargin < 4 || isempty(PopulationSize)
        PopulationSize = numel(Population);
    end

    if PopulationSize > numel(Population)
        error('PopulationSize is larger than actual population size.');
    end

    for i = 1:PopulationSize

        if ~isfield(Population(i), 'Obj') || numel(Population(i).Obj) < 2
            error('Population(%d) does not contain a valid Obj field.', i);
        end

        obj = Population(i).Obj;

        % Skip invalid objective values
        if any(isnan(obj)) || any(isinf(obj))
            continue;
        end

        % =========================
        % Objective 1: Conductance
        % =========================
        if obj(1) < IdealPoint(1)

            IdealPoint(1) = obj(1);

            IndivPoint(1).Chromosome = Population(i).Chromosome;
            IndivPoint(1).CmplxID = Population(i).CommsNodes;
            IndivPoint(1).Ideal = IdealPoint(1);
            IndivPoint(1).Obj = obj;

            if isfield(Population(i), 'Kv')
                IndivPoint(1).Kv = Population(i).Kv;
            end

            if isfield(Population(i), 'avgGO')
                IndivPoint(1).avgGO = Population(i).avgGO;
            end
        end

        % =========================
        % Objective 2: GO loss
        % =========================
        if obj(2) < IdealPoint(2)

            IdealPoint(2) = obj(2);

            IndivPoint(2).Chromosome = Population(i).Chromosome;
            IndivPoint(2).CmplxID = Population(i).CommsNodes;
            IndivPoint(2).Ideal = IdealPoint(2);
            IndivPoint(2).Obj = obj;

            if isfield(Population(i), 'Kv')
                IndivPoint(2).Kv = Population(i).Kv;
            end

            if isfield(Population(i), 'avgGO')
                IndivPoint(2).avgGO = Population(i).avgGO;
            end
        end
    end
end
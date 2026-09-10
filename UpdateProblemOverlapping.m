function [Individual] = UpdateProblemOverlapping(Individual, Child, ...
    SubProblems, Params, IndivPoint, ObjectiveDimension)
% UpdateProblemOverlapping
% MOEA/D neighborhood update using Tchebycheff scalarization.
%
% Minimization problem:
%   Obj(1): Conductance
%   Obj(2): GO coherence loss

    if ObjectiveDimension ~= 2
        error('UpdateProblemOverlapping currently supports 2 objectives only.');
    end

    for ProblemCounter = 1:Params.PopSize

        for NeighbourCounter = 1:Params.Neighbour

            NeighbourIndex = SubProblems(ProblemCounter).Neighbour(NeighbourCounter);

            IndividualOfNeighbour = Individual(NeighbourIndex);
            WeightOfNeighbour = SubProblems(NeighbourIndex).Weight;

            F1 = ScalarFunction( ...
                IndividualOfNeighbour, ...
                WeightOfNeighbour, ...
                IndivPoint, ...
                ObjectiveDimension);

            F2 = ScalarFunction( ...
                Child(ProblemCounter), ...
                WeightOfNeighbour, ...
                IndivPoint, ...
                ObjectiveDimension);

            if F2 < F1

                Individual(NeighbourIndex).Chromosome = Child(ProblemCounter).Chromosome;
                Individual(NeighbourIndex).CommsNodes = Child(ProblemCounter).CommsNodes;
                Individual(NeighbourIndex).Obj = Child(ProblemCounter).Obj;

                Individual(NeighbourIndex).Obj1_Conductance = Child(ProblemCounter).Obj(1);
                Individual(NeighbourIndex).Obj2_GOCoherenceLoss = Child(ProblemCounter).Obj(2);

                if isfield(Child(ProblemCounter), 'Kv')
                    Individual(NeighbourIndex).Kv = Child(ProblemCounter).Kv;
                end

                if isfield(Child(ProblemCounter), 'avgGO')
                    Individual(NeighbourIndex).avgGO = Child(ProblemCounter).avgGO;
                end
            end
        end
    end
end
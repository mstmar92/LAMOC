function [Population] = CreatePopulationOverlapping(proteinPairs, PopulationSize)
% CreatePopulationOverlapping
% Creates initial edge-based chromosomes for overlapping complex detection.
%
% proteinPairs: E x 2 matrix of edge endpoints.
% Population(i).Chromosome is a vector of length E.
% Each gene stores an index of a randomly selected adjacent edge.
%
% Note:
% Randomness is controlled externally by rng(...) in Main_TNSE.

    % Ensure valid indexing type
    proteinPairs = double(proteinPairs);

    % ---------------------------------
    % Basic sizes
    % ---------------------------------
    numEdges = size(proteinPairs, 1);
    numNodes = max(proteinPairs(:));

    % ---------------------------------
    % Build incident edge list for each node
    % incidentEdges{node} = list of edge indices touching that node
    % ---------------------------------
    incidentEdges = cell(numNodes, 1);

    for e = 1:numEdges
        u = proteinPairs(e, 1);
        v = proteinPairs(e, 2);

        incidentEdges{u}(end+1) = e;
        incidentEdges{v}(end+1) = e;
    end

    % ---------------------------------
    % Precompute candidate adjacent edges for each edge
    % ---------------------------------
    adjacentEdgeList = cell(numEdges, 1);

    for e = 1:numEdges
        u = proteinPairs(e, 1);
        v = proteinPairs(e, 2);

        candidates = [incidentEdges{u}, incidentEdges{v}];
        candidates = unique(candidates);
        candidates(candidates == e) = [];

        adjacentEdgeList{e} = candidates;
    end

    % ---------------------------------
    % Preallocate population
    % ---------------------------------
    Population = struct('Chromosome', cell(1, PopulationSize));

    for p = 1:PopulationSize

        chromosome = zeros(1, numEdges);

        for e = 1:numEdges

            candidates = adjacentEdgeList{e};

            if ~isempty(candidates)
                chromosome(e) = candidates(randi(numel(candidates)));
            else
                chromosome(e) = 0;
            end
        end

        Population(p).Chromosome = chromosome;
    end
end
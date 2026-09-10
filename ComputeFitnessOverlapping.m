function [Population] = ComputeFitnessOverlapping(A_ppi, A_go, N, ...
    IndicesInteractionProtein, NumInteractionProtein, MaxNumInteractionProtein, ...
    Population, PopulationSize)
% ComputeFitnessOverlapping
% Computes the two-objective fitness for overlapping protein complexes.
%
% Objectives:
%   Obj1 = average PPI conductance       (minimize)
%   Obj2 = 1 - average GO coherence      (minimize)
%
% Notes:
%   IndicesInteractionProtein, NumInteractionProtein, and
%   MaxNumInteractionProtein are kept in the signature for compatibility
%   with the original LAMOC pipeline, but are not used here.

    %#ok<*INUSD>

    if nargin < 8 || isempty(PopulationSize)
        PopulationSize = numel(Population);
    end

    if PopulationSize > numel(Population)
        error('PopulationSize is larger than the actual number of individuals.');
    end

    % ---------------------------------
    % Precompute static PPI network info
    % ---------------------------------
    A_ppi = sparse(double(A_ppi));

    if size(A_ppi,1) ~= N || size(A_ppi,2) ~= N
        error('A_ppi size does not match N.');
    end

    A_ppi(1:N+1:end) = 0;

    deg = full(sum(A_ppi, 2));
    totalVol = sum(deg);

    % ---------------------------------
    % Prepare GO layer
    % ---------------------------------
    if isempty(A_go) || nnz(A_go) == 0

        A_go = sparse(N, N);
        useGO = false;

    else

        A_go = sparse(double(A_go));

        if size(A_go,1) ~= N || size(A_go,2) ~= N
            error('A_go size does not match N.');
        end

        A_go(A_go < 0) = 0;
        A_go(1:N+1:end) = 0;

        useGO = true;
    end

    % ---------------------------------
    % Initialize output fields
    % ---------------------------------
    for i = 1:PopulationSize
        Population(i).Obj = [];
        Population(i).Obj1_Conductance = 0;
        Population(i).Obj2_GOCoherenceLoss = 0;
        Population(i).Kv = 0;
        Population(i).avgGO = 0;
    end

    % ---------------------------------
    % Compute fitness per individual
    % ---------------------------------
    for i = 1:PopulationSize

        if ~isfield(Population(i), 'CommsNodes')
            error('Population(%d) does not contain CommsNodes. Decode before fitness.', i);
        end

        Population(i) = ComputeOurFitnessCollectionOverlappingFast( ...
            A_ppi, A_go, deg, totalVol, useGO, Population(i));
    end
end
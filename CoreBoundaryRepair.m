function Child = CoreBoundaryRepair(Child, A_ppi, A_go, ...
                                    IndicesInteractionProtein, NumInteractionProtein, params)
% CoreBoundaryRepair
% Faster true core-boundary repair for overlapping protein complexes.
%
% Steps:
%   1) Score current members by topology, triangle support, and GO.
%   2) Select core nodes.
%   3) Remove weak boundary nodes.
%   4) Add only a limited number of promising boundary candidates.
%
% This version is optimized by limiting the candidate pool.

    % =========================
    % Default parameters
    % =========================
    if nargin < 6 || isempty(params)
        params = struct();
    end

    if ~isfield(params,'MinSize'), params.MinSize = 3; end
    if ~isfield(params,'CoreThreshold'), params.CoreThreshold = 0.60; end
    if ~isfield(params,'RemoveThreshold'), params.RemoveThreshold = 0.25; end

    if ~isfield(params,'AlphaTopo'), params.AlphaTopo = 0.50; end
    if ~isfield(params,'BetaTri'),   params.BetaTri   = 0.30; end
    if ~isfield(params,'GammaGO'),   params.GammaGO   = 0.20; end

    if ~isfield(params,'AddThreshold'), params.AddThreshold = 0.50; end
    if ~isfield(params,'GOAddThreshold'), params.GOAddThreshold = 0.10; end
    if ~isfield(params,'MinInternalNeighbors'), params.MinInternalNeighbors = 2; end
    if ~isfield(params,'MaxBoundaryAdds'), params.MaxBoundaryAdds = 2; end
    if ~isfield(params,'MaxCandidatePool'), params.MaxCandidatePool = 20; end

    A_ppi = sparse(double(A_ppi));
    A_go  = sparse(double(A_go));

    N = size(A_ppi,1);

    if ~isfield(Child, 'CommsNodes')
        error('Child must contain CommsNodes before CoreBoundaryRepair.');
    end

    Comms = Child.CommsNodes;
    K = numel(Comms);

    if K == 0
        return;
    end

    % =========================
    % Repair each complex
    % =========================
    for c = 1:K

        members = unique(double(Comms{c}(:)))';
        members = members(members >= 1 & members <= N);

        if numel(members) < params.MinSize
            Comms{c} = members;
            continue;
        end

        % ----------------------------------------------------
        % 1) Score current members
        % ----------------------------------------------------
        [affinity, topoScore, triScore, goScore] = LocalMemberAffinity( ...
            members, A_ppi, A_go, ...
            IndicesInteractionProtein, NumInteractionProtein, ...
            params, N);

        % ----------------------------------------------------
        % 2) Core selection
        % ----------------------------------------------------
        coreNodes = members(affinity >= params.CoreThreshold);

        if numel(coreNodes) < params.MinSize
            [~, order] = sort(affinity, 'descend');
            coreNodes = members(order(1:min(params.MinSize, numel(members))));
        end

        boundaryNodes = setdiff(members, coreNodes, 'stable');

        % ----------------------------------------------------
        % 3) Boundary pruning
        % ----------------------------------------------------
        keptBoundary = [];

        for b = 1:numel(boundaryNodes)
            u = boundaryNodes(b);
            idx = find(members == u, 1);

            if ~isempty(idx) && affinity(idx) >= params.RemoveThreshold
                keptBoundary(end+1) = u; %#ok<AGROW>
            end
        end

        repairedMembers = unique([coreNodes, keptBoundary], 'stable');

        % Preserve minimum size using best original members
        if numel(repairedMembers) < params.MinSize
            [~, order] = sort(affinity, 'descend');

            for oi = 1:numel(order)
                u = members(order(oi));

                if ~ismember(u, repairedMembers)
                    repairedMembers(end+1) = u; %#ok<AGROW>
                end

                if numel(repairedMembers) >= params.MinSize
                    break;
                end
            end
        end

        repairedMembers = unique(repairedMembers, 'stable');

        % ----------------------------------------------------
        % 4) Fast boundary expansion
        % ----------------------------------------------------
        if numel(coreNodes) >= params.MinSize && params.MaxBoundaryAdds > 0

            candidatePool = [];

            % collect neighbors of core nodes only
            for ii = 1:numel(coreNodes)

                u = coreNodes(ii);

                neigh = IndicesInteractionProtein(u, 1:NumInteractionProtein(u));
                neigh = double(neigh);
                neigh = neigh(neigh >= 1 & neigh <= N);

                candidatePool = [candidatePool, neigh]; %#ok<AGROW>
            end

            candidatePool = unique(candidatePool, 'stable');
            candidatePool = setdiff(candidatePool, repairedMembers, 'stable');

            if ~isempty(candidatePool)

                mask = false(N,1);
                mask(repairedMembers) = true;

                % quick topological pre-score: number of links to repairedMembers
                ninQuick = zeros(numel(candidatePool),1);

                for q = 1:numel(candidatePool)

                    u = candidatePool(q);

                    neigh = IndicesInteractionProtein(u, 1:NumInteractionProtein(u));
                    neigh = double(neigh);
                    neigh = neigh(neigh >= 1 & neigh <= N);

                    if isempty(neigh)
                        ninQuick(q) = 0;
                    else
                        ninQuick(q) = sum(mask(neigh));
                    end
                end

                % keep only candidates with enough internal neighbors
                validQuick = ninQuick >= params.MinInternalNeighbors;

                candidatePool = candidatePool(validQuick);
                ninQuick = ninQuick(validQuick);

                if ~isempty(candidatePool)

                    % limit candidate pool to top candidates
                    [~, ordQuick] = sort(ninQuick, 'descend');

                    keepCount = min(params.MaxCandidatePool, numel(ordQuick));

                    candidatePool = candidatePool(ordQuick(1:keepCount));

                    candScores = zeros(numel(candidatePool),1);
                    candGO = zeros(numel(candidatePool),1);
                    candNin = zeros(numel(candidatePool),1);

                    for q = 1:numel(candidatePool)

                        u = candidatePool(q);

                        [scoreU, goU, ninU] = LocalCandidateAffinity( ...
                            u, repairedMembers, A_ppi, A_go, ...
                            IndicesInteractionProtein, NumInteractionProtein, ...
                            params, N);

                        candScores(q) = scoreU;
                        candGO(q) = goU;
                        candNin(q) = ninU;
                    end

                    validAdd = candScores >= params.AddThreshold & ...
                               candGO >= params.GOAddThreshold & ...
                               candNin >= params.MinInternalNeighbors;

                    validCandidates = candidatePool(validAdd);
                    validScores = candScores(validAdd);

                    if ~isempty(validCandidates)

                        [~, ord] = sort(validScores, 'descend');

                        maxAdds = min(params.MaxBoundaryAdds, numel(ord));
                        addNodes = validCandidates(ord(1:maxAdds));

                        repairedMembers = unique([repairedMembers, addNodes], 'stable');
                    end
                end
            end
        end

        repairedMembers = unique(double(repairedMembers(:)))';
        repairedMembers = repairedMembers(repairedMembers >= 1 & repairedMembers <= N);

        Comms{c} = repairedMembers;
    end

    Child.CommsNodes = Comms;
end


% ========================================================================
% Local member affinity
% ========================================================================
function [affinity, topoScore, triScore, goScore] = LocalMemberAffinity( ...
    members, A_ppi, A_go, ...
    IndicesInteractionProtein, NumInteractionProtein, ...
    params, N)

    members = unique(double(members(:)))';

    mask = false(N,1);
    mask(members) = true;

    m = numel(members);

    affinity = zeros(1,m);
    topoScore = zeros(1,m);
    triScore = zeros(1,m);
    goScore = zeros(1,m);

    complexSize = numel(members);

    for ii = 1:m

        u = members(ii);

        neigh = IndicesInteractionProtein(u, 1:NumInteractionProtein(u));
        neigh = double(neigh);
        neigh = neigh(neigh >= 1 & neigh <= N);

        if isempty(neigh)
            continue;
        end

        inMask = mask(neigh);

        kin = sum(inMask);
        kout = numel(neigh) - kin;

        topoScore(ii) = kin / (kin + kout + 1e-6);

        % Faster triangle support approximation:
        % For internal neighbors v, count how many of v's neighbors are in the complex.
        triSupport = 0;

        internalNeigh = neigh(inMask);

        for t = 1:numel(internalNeigh)

            v = internalNeigh(t);

            vn = IndicesInteractionProtein(v, 1:NumInteractionProtein(v));
            vn = double(vn);
            vn = vn(vn >= 1 & vn <= N);

            triSupport = triSupport + sum(mask(vn));
        end

        triScore(ii) = triSupport / (max(1, kin) * max(1, complexSize));
        triScore(ii) = min(triScore(ii), 1);

        others = members(members ~= u);

        if isempty(others)
            goScore(ii) = 0;
        else
            goScore(ii) = full(mean(A_go(u, others)));
        end

        goScore(ii) = min(max(goScore(ii), 0), 1);

        affinity(ii) = params.AlphaTopo * topoScore(ii) + ...
                       params.BetaTri   * triScore(ii) + ...
                       params.GammaGO   * goScore(ii);

        affinity(ii) = min(max(affinity(ii), 0), 1);
    end
end


% ========================================================================
% Local candidate affinity
% ========================================================================
function [affinity, goScore, nin] = LocalCandidateAffinity( ...
    u, complexMembers, A_ppi, A_go, ...
    IndicesInteractionProtein, NumInteractionProtein, ...
    params, N)

    complexMembers = unique(double(complexMembers(:)))';

    mask = false(N,1);
    mask(complexMembers) = true;

    neigh = IndicesInteractionProtein(u, 1:NumInteractionProtein(u));
    neigh = double(neigh);
    neigh = neigh(neigh >= 1 & neigh <= N);

    if isempty(neigh)
        affinity = 0;
        goScore = 0;
        nin = 0;
        return;
    end

    inMask = mask(neigh);

    nin = sum(inMask);
    nout = numel(neigh) - nin;

    topoScore = nin / (nin + nout + 1e-6);

    % Fast triangle support
    triSupport = 0;

    internalNeigh = neigh(inMask);

    for t = 1:numel(internalNeigh)

        v = internalNeigh(t);

        vn = IndicesInteractionProtein(v, 1:NumInteractionProtein(v));
        vn = double(vn);
        vn = vn(vn >= 1 & vn <= N);

        triSupport = triSupport + sum(mask(vn));
    end

    triScore = triSupport / (max(1, nin) * max(1, numel(complexMembers)));
    triScore = min(triScore, 1);

    if isempty(complexMembers)
        goScore = 0;
    else
        goScore = full(mean(A_go(u, complexMembers)));
    end

    goScore = min(max(goScore, 0), 1);

    affinity = params.AlphaTopo * topoScore + ...
               params.BetaTri   * triScore + ...
               params.GammaGO   * goScore;

    affinity = min(max(affinity, 0), 1);
end
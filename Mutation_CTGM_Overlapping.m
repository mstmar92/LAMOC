function Child = Mutation_CTGM_Overlapping(Child, A_ppi, A_go, ...
                                           IndicesInteractionProtein, NumInteractionProtein, ...
                                           Pm, params)
% Mutation_CTGM_Overlapping
% Heuristic-guided overlapping mutation using:
%   - PPI topology
%   - GO semantic similarity
%
% Important:
% This mutation operates on Child.CommsNodes, not on Child.Chromosome.
% Therefore, it must be called AFTER decoding the chromosome.

    if nargin < 7 || isempty(params)
        params = struct();
    end

    if ~isfield(params,'MinSize'), params.MinSize = 3; end
    if ~isfield(params,'pAddTriangle'), params.pAddTriangle = 0.6; end
    if ~isfield(params,'pPrune'), params.pPrune = 0.4; end
    if ~isfield(params,'MaxComplexAddsPerNode'), params.MaxComplexAddsPerNode = 2; end
    if ~isfield(params,'MaxMovesPerNode'), params.MaxMovesPerNode = 1; end
    if ~isfield(params,'GOScoreWeight'), params.GOScoreWeight = 1.0; end
    if ~isfield(params,'GOAddThreshold'), params.GOAddThreshold = 0.15; end
    if ~isfield(params,'GOPruneThreshold'), params.GOPruneThreshold = 0.10; end

    A_ppi = sparse(double(A_ppi));
    A_go  = sparse(double(A_go));

    N = size(A_ppi,1);

    if ~isfield(Child, 'CommsNodes')
        error('Child must contain CommsNodes before Mutation_CTGM_Overlapping.');
    end

    Comms = Child.CommsNodes;
    K = numel(Comms);

    if K == 0
        return;
    end

    % =========================
    % Build membership map
    % =========================
    proteinComplexMap = cell(N,1);

    for c = 1:K

        members = unique(double(Comms{c}(:)))';
        members = members(members >= 1 & members <= N);

        Comms{c} = members;

        for ii = 1:numel(members)
            u = members(ii);
            proteinComplexMap{u}(end+1) = c;
        end
    end

    % =========================
    % Build complex masks
    % =========================
    Cmask = false(N, K);

    for c = 1:K
        if ~isempty(Comms{c})
            Cmask(Comms{c}, c) = true;
        end
    end

    % =========================
    % Mutation over nodes
    % =========================
    for u = 1:N

        if NumInteractionProtein(u) <= 0 || rand > Pm
            continue;
        end

        currCs = proteinComplexMap{u};

        if isempty(currCs)
            continue;
        end

        neigh = IndicesInteractionProtein(u, 1:NumInteractionProtein(u));
        neigh = double(neigh);
        neigh = neigh(neigh >= 1 & neigh <= N);

        if isempty(neigh)
            continue;
        end

        % ==========================================================
        % A) Boundary-aware add/move
        % ==========================================================
        moveCount = 0;

        for idxC = 1:numel(currCs)

            if moveCount >= params.MaxMovesPerNode
                break;
            end

            c = currCs(idxC);

            if c < 1 || c > K
                continue;
            end

            membersMask = Cmask(:,c);

            kin = 0;
            kout = 0;

            for t = 1:numel(neigh)

                v = neigh(t);

                if membersMask(v)
                    kin = kin + A_ppi(u,v);
                else
                    kout = kout + A_ppi(u,v);
                end
            end

            if kin <= kout

                bestC = [];
                bestScore = -Inf;

                cand = [];

                for t = 1:numel(neigh)
                    v = neigh(t);
                    cand = [cand, proteinComplexMap{v}]; %#ok<AGROW>
                end

                cand = unique(cand);
                cand(cand == c) = [];
                cand = cand(cand >= 1 & cand <= K);

                for cc = cand

                    if Cmask(u,cc)
                        continue;
                    end

                    mask2 = Cmask(:,cc);
                    members2 = find(mask2);

                    kin2 = 0;
                    kout2 = 0;

                    for t = 1:numel(neigh)

                        v = neigh(t);

                        if mask2(v)
                            kin2 = kin2 + A_ppi(u,v);
                        else
                            kout2 = kout2 + A_ppi(u,v);
                        end
                    end

                    triGain = 0;

                    for t = 1:numel(neigh)

                        v = neigh(t);

                        if ~mask2(v) || A_ppi(u,v) == 0
                            continue;
                        end

                        vn = IndicesInteractionProtein(v, 1:NumInteractionProtein(v));
                        vn = double(vn);
                        vn = vn(vn >= 1 & vn <= N);

                        triGain = triGain + sum(mask2(vn));
                    end

                    if isempty(members2)
                        goGain = 0;
                    else
                        goGain = full(mean(A_go(u, members2)));
                    end

                    score = full(kin2 - kout2) + ...
                            0.1 * triGain + ...
                            params.GOScoreWeight * goGain;

                    if score > bestScore
                        bestScore = score;
                        bestC = cc;
                    end
                end

                if ~isempty(bestC)

                    Comms{bestC}(end+1) = u;
                    Comms{bestC} = unique(Comms{bestC});

                    Cmask(u,bestC) = true;

                    proteinComplexMap{u}(end+1) = bestC;
                    proteinComplexMap{u} = unique(proteinComplexMap{u});

                    moveCount = moveCount + 1;
                end
            end
        end

        % ==========================================================
        % B) Triangle-closure add with GO support
        % ==========================================================
        if rand <= params.pAddTriangle

            adds = 0;
            cand = [];

            for t = 1:numel(neigh)
                v = neigh(t);
                cand = [cand, proteinComplexMap{v}]; %#ok<AGROW>
            end

            cand = unique(cand);
            cand = cand(cand >= 1 & cand <= K);

            for cc = cand

                if adds >= params.MaxComplexAddsPerNode
                    break;
                end

                if Cmask(u,cc)
                    continue;
                end

                mask = Cmask(:,cc);
                members_cc = find(mask);

                nin = 0;

                for t = 1:numel(neigh)

                    v = neigh(t);

                    if mask(v) && A_ppi(u,v) ~= 0
                        nin = nin + 1;

                        if nin >= 2
                            break;
                        end
                    end
                end

                if isempty(members_cc)
                    goGain = 0;
                else
                    goGain = full(mean(A_go(u, members_cc)));
                end

                if nin >= 2 && goGain >= params.GOAddThreshold

                    Comms{cc}(end+1) = u;
                    Comms{cc} = unique(Comms{cc});

                    Cmask(u,cc) = true;

                    proteinComplexMap{u}(end+1) = cc;
                    proteinComplexMap{u} = unique(proteinComplexMap{u});

                    adds = adds + 1;
                end
            end
        end

        % ==========================================================
        % C) Prune weak boundary nodes
        % ==========================================================
        if rand <= params.pPrune

            currCs2 = unique(proteinComplexMap{u});
            currCs2 = currCs2(currCs2 >= 1 & currCs2 <= K);

            for idxC = 1:numel(currCs2)

                c = currCs2(idxC);
                members = unique(double(Comms{c}(:)))';

                if numel(members) <= params.MinSize
                    continue;
                end

                mask = Cmask(:,c);

                kin = 0;
                kout = 0;
                triSupport = 0;

                for t = 1:numel(neigh)

                    v = neigh(t);

                    if mask(v)

                        kin = kin + A_ppi(u,v);

                        vn = IndicesInteractionProtein(v, 1:NumInteractionProtein(v));
                        vn = double(vn);
                        vn = vn(vn >= 1 & vn <= N);

                        triSupport = triSupport + sum(mask(vn));

                    else
                        kout = kout + A_ppi(u,v);
                    end
                end

                members_wo_u = members(members ~= u);

                if isempty(members_wo_u)
                    goSupport = 0;
                else
                    goSupport = full(mean(A_go(u, members_wo_u)));
                end

                if full(kin <= kout) && ...
                   triSupport <= 1 && ...
                   goSupport < params.GOPruneThreshold

                    Comms{c} = members(members ~= u);
                    Cmask(u,c) = false;

                    proteinComplexMap{u} = proteinComplexMap{u}(proteinComplexMap{u} ~= c);
                end
            end
        end
    end

    % =========================
    % Final clean
    % =========================
    for c = 1:K
        Comms{c} = unique(double(Comms{c}(:)))';
        Comms{c} = Comms{c}(Comms{c} >= 1 & Comms{c} <= N);
    end

    Child.CommsNodes = Comms;
end
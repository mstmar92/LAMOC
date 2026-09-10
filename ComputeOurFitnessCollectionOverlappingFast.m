function [Population] = ComputeOurFitnessCollectionOverlappingFast( ...
    A_ppi, A_go, deg, totalVol, useGO, Population)
% ComputeOurFitnessCollectionOverlappingFast
% Computes fitness for one individual.
%
% Obj1 = average conductance over predicted complexes.
% Obj2 = 1 - average GO coherence.
%
% Both objectives are minimized.

    MinSize = 3;
    epsVal  = 1e-12;

    N = size(A_ppi, 1);

    if ~isfield(Population, 'CommsNodes')
        error('Population does not contain CommsNodes.');
    end

    CmplxID = Population.CommsNodes;
    K = numel(CmplxID);

    condSum = 0;
    goSum   = 0;
    Kv = 0;

    for c = 1:K

        members = CmplxID{c};

        if isempty(members)
            continue;
        end

        members = unique(double(members(:)))';

        % Safety: keep only valid protein indices
        members = members(members >= 1 & members <= N);

        s = numel(members);

        if s < MinSize
            continue;
        end

        Kv = Kv + 1;

        % =========================
        % Obj1: Conductance on PPI
        % =========================
        A_sub_ppi = A_ppi(members, members);

        m_in = nnz(triu(A_sub_ppi, 1));

        volC = sum(deg(members));
        cutC = volC - 2 * m_in;

        denom = min(volC, totalVol - volC);

        if denom <= epsVal
            phi = 1;
        else
            phi = cutC / (denom + epsVal);
        end

        % Numerical safety
        if phi < 0
            phi = 0;
        elseif phi > 1
            phi = 1;
        end

        condSum = condSum + phi;

        % =========================
        % Obj2: GO coherence
        % =========================
        if useGO

            A_sub_go = A_go(members, members);

            pairCount = s * (s - 1) / 2;

            if pairCount <= 0
                goCoh = 0;
            else
                goCoh = full(sum(sum(triu(A_sub_go, 1)))) / (pairCount + epsVal);
            end

            % Wang similarity should be in [0,1]
            if goCoh < 0
                goCoh = 0;
            elseif goCoh > 1
                goCoh = 1;
            end

        else
            goCoh = 0;
        end

        goSum = goSum + goCoh;
    end

    % =========================
    % Aggregate objectives
    % =========================
    if Kv == 0

        Obj1 = 1;
        avgGO = 0;
        Obj2 = 1;

    else

        Obj1 = condSum / Kv;
        avgGO = goSum / Kv;
        Obj2 = 1 - avgGO;

        % Safety bounds
        Obj1 = min(max(Obj1, 0), 1);
        avgGO = min(max(avgGO, 0), 1);
        Obj2 = min(max(Obj2, 0), 1);
    end

    Population.Obj = [Obj1, Obj2];
    Population.Obj1_Conductance = Obj1;
    Population.Obj2_GOCoherenceLoss = Obj2;
    Population.Kv = Kv;
    Population.avgGO = avgGO;
end
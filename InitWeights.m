function SubProblem = InitWeights(PopSize, Neighbour, ObjectiveDimension)
% InitWeights
% Initializes MOEA/D decomposition weight vectors and neighborhood structure.
%
% For ObjectiveDimension = 2:
%   Generates PopSize evenly spaced weights from [0,1] to [1,0].
%
% Output:
%   SubProblem(i).Weight
%   SubProblem(i).Neighbour
%   SubProblem(i).Optimal
%   SubProblem(i).OptPoint
%   SubProblem(i).Curpoint

    if nargin < 3
        error('InitWeights requires PopSize, Neighbour, and ObjectiveDimension.');
    end

    if PopSize < 2
        error('PopSize must be at least 2.');
    end

    Neighbour = min(Neighbour, PopSize);

    % =========================
    % Generate weight vectors
    % =========================
    if ObjectiveDimension == 2

        SubProblem = repmat( ...
            struct('Weight', [], 'Neighbour', [], 'Optimal', Inf, ...
                   'OptPoint', [], 'Curpoint', []), ...
            1, PopSize);

        for i = 1:PopSize
            w1 = (i - 1) / (PopSize - 1);
            w2 = 1 - w1;

            SubProblem(i).Weight = [w1, w2];
        end

    elseif ObjectiveDimension == 3

        temp = struct('Weight', {}, 'Neighbour', {}, 'Optimal', {}, ...
                      'OptPoint', {}, 'Curpoint', {});

        counter = 0;

        for i = 0:PopSize
            for j = 0:(PopSize - i)
                k = PopSize - i - j;

                counter = counter + 1;

                temp(counter).Weight = [i, j, k] / PopSize;
                temp(counter).Neighbour = [];
                temp(counter).Optimal = Inf;
                temp(counter).OptPoint = [];
                temp(counter).Curpoint = [];
            end
        end

        SubProblem = temp;

    else
        error('Only ObjectiveDimension = 2 or 3 is supported.');
    end

    % =========================
    % Build neighborhood
    % =========================
    Length = numel(SubProblem);
    Neighbour = min(Neighbour, Length);

    WeightMatrix = zeros(Length, ObjectiveDimension);

    for i = 1:Length
        WeightMatrix(i,:) = SubProblem(i).Weight;
    end

    DistanceMatrix = zeros(Length, Length);

    for i = 1:Length
        for j = i+1:Length
            d = WeightMatrix(i,:) - WeightMatrix(j,:);
            DistanceMatrix(i,j) = d * d';
            DistanceMatrix(j,i) = DistanceMatrix(i,j);
        end
    end

    for i = 1:Length
        [~, sindex] = sort(DistanceMatrix(i,:));
        SubProblem(i).Neighbour = sindex(1:Neighbour)';
    end
end
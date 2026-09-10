function [Results] = MOEAD_Overlapping_TNSE(A_ppi, A_go, N, ...
    IndicesInteractionProtein, NumInteractionProtein, MaxNumInteractionProtein, ...
    Population, PopulationSize, ...
    Heuristic, Pm, ...
    ChromHeuristic, ProteinLabelPairs, ...
    UseGOinFitness, UseHeuristicMutation, UseRepair, ...
    Generations)

%-------------------------------------------------------------------------%
%                    MOEA/D Overlapping - TNSE Version                     %
%-------------------------------------------------------------------------%
% Official lightweight version for TNSE-style experiments.
%
% Stores:
%   - FinalNearParetoOptimalSet
%   - FinalPopulation
%   - Final ideal point
%   - Runtime logs
%   - Convergence curves
%
% Does NOT store full Pareto set for every generation.
%-------------------------------------------------------------------------%

% =========================
% Default generations
% =========================
if nargin < 18 || isempty(Generations)
    Generations = 100;
end

%-------------------------------------------------------------------------%
%                        MOEA/D Parameter Settings                         %
%-------------------------------------------------------------------------%

ObjectiveDimension = 2;

IdealPoint = 1e7 * ones(1, ObjectiveDimension);
IndivPoint = [];

Params.PopSize = PopulationSize;
Params.Neighbour = min(5, PopulationSize);
Params.Generations = Generations;
Params.DMethod = 'ts';

SubProblem = InitWeights(Params.PopSize, Params.Neighbour, ObjectiveDimension);
Params.PopSize = length(SubProblem);

% =========================
% Mutation parameters
% =========================
Params.Mut.MinSize = 3;
Params.Mut.pAddTriangle = 0.6;
Params.Mut.pPrune = 0.4;
Params.Mut.MaxComplexAddsPerNode = 2;
Params.Mut.MaxMovesPerNode = 1;

Params.Mut.GOScoreWeight    = 1.0;
Params.Mut.GOAddThreshold   = 0.15;
Params.Mut.GOPruneThreshold = 0.10;

% =========================
% Repair parameters
% =========================
Params.Repair.MinSize = 3;

% Core-boundary classification
Params.Repair.CoreThreshold = 0.60;
Params.Repair.RemoveThreshold = 0.25;

% Affinity weights
Params.Repair.AlphaTopo = 0.50;
Params.Repair.BetaTri   = 0.30;
Params.Repair.GammaGO   = 0.20;

% Boundary expansion
Params.Repair.AddThreshold = 0.45;
Params.Repair.GOAddThreshold = 0.10;
Params.Repair.MinInternalNeighbors = 2;
Params.Repair.MaxBoundaryAdds = 5;

% =========================
% Prepare matrices
% =========================
A_ppi = sparse(double(A_ppi));

if UseGOinFitness || UseHeuristicMutation || UseRepair
    A_go = sparse(double(A_go));
else
    A_go = sparse(size(A_ppi,1), size(A_ppi,2));
end

% =========================
% GO matrix used in fitness
% =========================
if UseGOinFitness
    AgoFitness = A_go;
else
    AgoFitness = sparse(size(A_ppi,1), size(A_ppi,2));
end

%-------------------------------------------------------------------------%
% Runtime and convergence logs
%-------------------------------------------------------------------------%

GenerationTime = zeros(Params.Generations,1);

TimeLog.Selection        = zeros(Params.Generations,1);
TimeLog.Crossover        = zeros(Params.Generations,1);
TimeLog.Decode           = zeros(Params.Generations,1);
TimeLog.Mutation         = zeros(Params.Generations,1);
TimeLog.Repair           = zeros(Params.Generations,1);
TimeLog.Fitness          = zeros(Params.Generations,1);
TimeLog.UpdateReference  = zeros(Params.Generations,1);
TimeLog.UpdatePopulation = zeros(Params.Generations,1);
TimeLog.UpdatePareto     = zeros(Params.Generations,1);

Convergence.BestObj1     = zeros(Params.Generations,1);
Convergence.BestObj2     = zeros(Params.Generations,1);
Convergence.MeanObj1     = zeros(Params.Generations,1);
Convergence.MeanObj2     = zeros(Params.Generations,1);
Convergence.ParetoSize   = zeros(Params.Generations,1);
Convergence.MeanKv       = zeros(Params.Generations,1);
Convergence.MeanAvgGO    = zeros(Params.Generations,1);

%-------------------------------------------------------------------------%
%                              Initialization                              %
%-------------------------------------------------------------------------%

t_init = tic;

[IdealPoint, IndivPoint] = UpdateReferenceOverlapping( ...
    IdealPoint, IndivPoint, Population, Params.PopSize);

NearParetoOptimalSet = [];

NearParetoOptimalSet = UpdateNearParetoOptimalSetOverlapping( ...
    NearParetoOptimalSet, Population, Params.PopSize, ObjectiveDimension);

GenerationTime(1) = toc(t_init);

% Initial convergence values
[Convergence] = UpdateConvergenceLog( ...
    Convergence, 1, Population, Params.PopSize, NearParetoOptimalSet);

fprintf('Initial reference + Pareto update: %.2f sec\n', GenerationTime(1));

%-------------------------------------------------------------------------%
%                              Evolution Loop                              %
%-------------------------------------------------------------------------%

for GenerationCounter = 2:Params.Generations

    printThisGeneration = mod(GenerationCounter,10) == 0 || ...
                          GenerationCounter == 2 || ...
                          GenerationCounter == Params.Generations;

    if printThisGeneration
        fprintf('\nGeneration %d / %d\n', GenerationCounter, Params.Generations);
    end

    t_gen = tic;

    Child = Population;

    % =========================
    % Selection
    % =========================
    t1 = tic;

    Parent1 = repmat(Population(1), 1, Params.PopSize);
    Parent2 = repmat(Population(1), 1, Params.PopSize);

    for ProblemCounter = 1:Params.PopSize

        rand1 = randi(Params.Neighbour);
        rand2 = randi(Params.Neighbour);

        P1 = SubProblem(ProblemCounter).Neighbour(rand1);
        P2 = SubProblem(ProblemCounter).Neighbour(rand2);

        Parent1(ProblemCounter) = Population(P1);
        Parent2(ProblemCounter) = Population(P2);
    end

    TimeLog.Selection(GenerationCounter) = toc(t1);

    % =========================
    % Crossover
    % =========================
    t2 = tic;

    for i = 1:Params.PopSize
        Child(i) = Crossover(Parent1(i), Parent2(i), N, Child(i));
    end

    TimeLog.Crossover(GenerationCounter) = toc(t2);

    % =========================
    % Decode after crossover
    % =========================
    t3 = tic;

    Child = Individual2CmplxDecodingOverlapping( ...
        ProteinLabelPairs, Child, Params.PopSize);

    TimeLog.Decode(GenerationCounter) = toc(t3);

    % =========================
    % Heuristic Mutation
    % Must be after decoding because it modifies CommsNodes
    % =========================
    t4 = tic;

    if UseHeuristicMutation
        for i = 1:Params.PopSize
            Child(i) = Mutation_CTGM_Overlapping( ...
                Child(i), A_ppi, A_go, ...
                IndicesInteractionProtein, NumInteractionProtein, ...
                Pm, Params.Mut);
        end
    end

    TimeLog.Mutation(GenerationCounter) = toc(t4);

    % =========================
    % Core-Boundary Repair
    % =========================
    t5 = tic;

    if UseRepair
        for i = 1:Params.PopSize
            Child(i) = CoreBoundaryRepair( ...
                Child(i), A_ppi, A_go, ...
                IndicesInteractionProtein, NumInteractionProtein, Params.Repair);
        end
    end

    TimeLog.Repair(GenerationCounter) = toc(t5);

    % =========================
    % Fitness computation
    % =========================
    t6 = tic;

    Child = ComputeFitnessOverlapping( ...
        A_ppi, AgoFitness, N, ...
        IndicesInteractionProtein, NumInteractionProtein, MaxNumInteractionProtein, ...
        Child, Params.PopSize);

    TimeLog.Fitness(GenerationCounter) = toc(t6);

    % =========================
    % Update reference point
    % =========================
    t7 = tic;

    [IdealPoint, IndivPoint] = UpdateReferenceOverlapping( ...
        IdealPoint, IndivPoint, Child, Params.PopSize);

    TimeLog.UpdateReference(GenerationCounter) = toc(t7);

    % =========================
    % Update population
    % =========================
    t8 = tic;

    Population = UpdateProblemOverlapping( ...
        Population, Child, SubProblem, Params, IndivPoint, ObjectiveDimension);

    TimeLog.UpdatePopulation(GenerationCounter) = toc(t8);

    % =========================
    % Update Pareto set
    % =========================
    t9 = tic;

    NearParetoOptimalSet = UpdateNearParetoOptimalSetOverlapping( ...
        NearParetoOptimalSet, Child, Params.PopSize, ObjectiveDimension);

    TimeLog.UpdatePareto(GenerationCounter) = toc(t9);

    % =========================
    % Runtime and convergence
    % =========================
    GenerationTime(GenerationCounter) = toc(t_gen);

    Convergence = UpdateConvergenceLog( ...
        Convergence, GenerationCounter, Population, Params.PopSize, NearParetoOptimalSet);

    if printThisGeneration
        fprintf('  Selection        : %.2f sec\n', TimeLog.Selection(GenerationCounter));
        fprintf('  Crossover        : %.2f sec\n', TimeLog.Crossover(GenerationCounter));
        fprintf('  Decode           : %.2f sec\n', TimeLog.Decode(GenerationCounter));
        fprintf('  Mutation         : %.2f sec\n', TimeLog.Mutation(GenerationCounter));
        fprintf('  Repair           : %.2f sec\n', TimeLog.Repair(GenerationCounter));
        fprintf('  Fitness          : %.2f sec\n', TimeLog.Fitness(GenerationCounter));
        fprintf('  UpdateReference  : %.2f sec\n', TimeLog.UpdateReference(GenerationCounter));
        fprintf('  UpdatePopulation : %.2f sec\n', TimeLog.UpdatePopulation(GenerationCounter));
        fprintf('  UpdatePareto     : %.2f sec\n', TimeLog.UpdatePareto(GenerationCounter));
        fprintf('  Total generation : %.2f sec\n', GenerationTime(GenerationCounter));
        fprintf('  Pareto size      : %d\n', Convergence.ParetoSize(GenerationCounter));
        fprintf('  Best Obj         : [%.4f, %.4f]\n', ...
            Convergence.BestObj1(GenerationCounter), ...
            Convergence.BestObj2(GenerationCounter));
    end
end

%-------------------------------------------------------------------------%
%                              Final Results                               %
%-------------------------------------------------------------------------%

Results = struct();

Results.FinalNearParetoOptimalSet = NearParetoOptimalSet;
Results.FinalPopulation = Population;
Results.FinalIdealPoint = IdealPoint;
Results.FinalIndivPoint = IndivPoint;

Results.Convergence = Convergence;

Results.Meta.Params = Params;
Results.Meta.ObjectiveDimension = ObjectiveDimension;

Results.Meta.UseGOinFitness = UseGOinFitness;
Results.Meta.UseHeuristicMutation = UseHeuristicMutation;
Results.Meta.UseRepair = UseRepair;

Results.Meta.GenerationTime = GenerationTime;
Results.Meta.TimeLog = TimeLog;

Results.Meta.TotalMOEADTime = sum(GenerationTime);

if Params.Generations > 1
    Results.Meta.MeanGenerationTime = mean(GenerationTime(2:end));
else
    Results.Meta.MeanGenerationTime = GenerationTime(1);
end

end


% ========================================================================
% Helper: Update convergence log
% ========================================================================
function Convergence = UpdateConvergenceLog( ...
    Convergence, gen, Population, PopSize, NearParetoOptimalSet)

    objMat = zeros(PopSize, 2);
    kvVals = zeros(PopSize, 1);
    goVals = zeros(PopSize, 1);

    for i = 1:PopSize

        if isfield(Population(i), 'Obj') && numel(Population(i).Obj) >= 2
            objMat(i,:) = Population(i).Obj(:)';
        else
            objMat(i,:) = [NaN, NaN];
        end

        if isfield(Population(i), 'Kv')
            kvVals(i) = Population(i).Kv;
        else
            kvVals(i) = NaN;
        end

        if isfield(Population(i), 'avgGO')
            goVals(i) = Population(i).avgGO;
        else
            goVals(i) = NaN;
        end
    end

    validObj = all(~isnan(objMat), 2);

    if any(validObj)
        Convergence.BestObj1(gen) = min(objMat(validObj,1));
        Convergence.BestObj2(gen) = min(objMat(validObj,2));
        Convergence.MeanObj1(gen) = mean(objMat(validObj,1));
        Convergence.MeanObj2(gen) = mean(objMat(validObj,2));
    else
        Convergence.BestObj1(gen) = NaN;
        Convergence.BestObj2(gen) = NaN;
        Convergence.MeanObj1(gen) = NaN;
        Convergence.MeanObj2(gen) = NaN;
    end

    Convergence.ParetoSize(gen) = numel(NearParetoOptimalSet);

    Convergence.MeanKv(gen) = mean(kvVals(~isnan(kvVals)));

    if isempty(goVals(~isnan(goVals)))
        Convergence.MeanAvgGO(gen) = NaN;
    else
        Convergence.MeanAvgGO(gen) = mean(goVals(~isnan(goVals)));
    end
end
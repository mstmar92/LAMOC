function Main_TNSE(NetworkNumber, ExpChoice, MaxRun, PopulationSize, Generations)
% Main_TNSE
% Official TNSE-ready experiment runner for LAMOC on FMvPCI yeast datasets.
%
% Usage:
%   Main_TNSE(6, 1, 30, 100, 100);  % Collins-FMvPCI TOPO
%   Main_TNSE(6, 2, 30, 100, 100);  % Collins-FMvPCI PPI_GO
%   Main_TNSE(6, 3, 30, 100, 100);  % Collins-FMvPCI FULL

    clc;

    % =========================
    % Default arguments
    % =========================
    if nargin < 3 || isempty(MaxRun)
        MaxRun = 30;
    end

    if nargin < 4 || isempty(PopulationSize)
        PopulationSize = 100;
    end

    if nargin < 5 || isempty(Generations)
        Generations = 100;
    end

    % =========================
    % Main parameters
    % =========================
    ChromHeuristic = 1;
    Heuristic = 1;
    Pm = 0.2;

    % No threshold for Wang semantic similarity
    T_go = 0;

    % =========================
    % Experiment type
    % =========================
    switch ExpChoice
        case 1
            ExperimentType = 'TOPO';
        case 2
            ExperimentType = 'PPI_GO';
        case 3
            ExperimentType = 'FULL';
        otherwise
            error('Invalid experiment type.');
    end

    % =========================
    % Dataset selection
    % =========================
    switch NetworkNumber

        case 1
            DatasetName = 'Yeast-D1';
            ProteinFile = 'DataSets\Protein\1-Protein-Yeast-D1-Files.mat';
            GOFile      = 'DataSets\GO\GO_Layer_YeastD1.mat';

        case 2
            DatasetName = 'Yeast-D2';
            ProteinFile = 'DataSets\Protein\2-Protein-Yeast-D2-Files.mat';
            GOFile      = 'DataSets\GO\GO_Layer_YeastD2.mat';

        case 3
            DatasetName = 'Collins-CYC';
            ProteinFile = 'DataSets\Protein\4-Collins.mat';
            GOFile      = 'DataSets\GO\GO_Layer_CollinsCYC.mat';

        case 4
            DatasetName = 'Collins-MIPS';
            ProteinFile = 'DataSets\Protein\5-Collins.mat';
            GOFile      = 'DataSets\GO\GO_Layer_CollinsMIPS.mat';

        case 5
            DatasetName = 'Human-BioGRID-Core';
            ProteinFile = 'DataSets\Protein\6-Human-BioGRID-Core-Files.mat';
            GOFile      = 'DataSets\GO\GO_Layer_HumanBioGRID_Core_FullResnik.mat';

        case 6
            DatasetName = 'Collins-FMvPCI';
            ProteinFile = 'DataSets\Protein\Collins_FMvPCI_Files.mat';
            GOFile      = 'DataSets\GO\GO_Layer_Collins_FMvPCI.mat';

        case 7
            DatasetName = 'Krogan-FMvPCI';
            ProteinFile = 'DataSets\Protein\Krogan_FMvPCI_Files.mat';
            GOFile      = 'DataSets\GO\GO_Layer_Krogan_FMvPCI.mat';

        case 8
            DatasetName = 'Gavin-FMvPCI';
            ProteinFile = 'DataSets\Protein\Gavin_FMvPCI_Files.mat';
            GOFile      = 'DataSets\GO\GO_Layer_Gavin_FMvPCI.mat';

        case 9
            DatasetName = 'DIP-Hsapi-FMvPCI';
            ProteinFile = 'DataSets\Protein\DIP-Hsapi_FMvPCI_Files_symbols.mat';
            GOFile      = 'DataSets\GO\GO_Layer_DIP-Hsapi_FMvPCI.mat';            

        otherwise
            error('Invalid network number.');
    end

    fprintf('\n=============================================\n');
    fprintf('LAMOC TNSE Experiment\n');
    fprintf('Dataset        : %s\n', DatasetName);
    fprintf('Experiment     : %s\n', ExperimentType);
    fprintf('MaxRun         : %d\n', MaxRun);
    fprintf('PopulationSize : %d\n', PopulationSize);
    fprintf('Generations    : %d\n', Generations);
    fprintf('=============================================\n');

    % =========================
    % Reproducibility metadata
    % =========================
    RunTimestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    MATLABVersion = version;

    GOType = 'Wang';
    GOOntologyFile = 'DataSets\GO\go-basic.obo';
    GOAnnotationFile = 'DataSets\GO\sgd.gaf';

    InitSeeds = zeros(MaxRun,1);
    RunSeeds  = zeros(MaxRun,1);

    for r = 1:MaxRun
        InitSeeds(r) = 1000 + r;
        RunSeeds(r)  = 100000 + r;
    end

    % =========================
    % Load PPI
    % =========================
    fprintf('\nLoading PPI file:\n%s\n', ProteinFile);

    if ~exist(ProteinFile, 'file')
        error('Protein file not found: %s', ProteinFile);
    end

    Sppi = load(ProteinFile);

    requiredVars = { ...
        'A', ...
        'N', ...
        'ProteinLabelPairs', ...
        'ProteinLabel', ...
        'NumInteractionProtein', ...
        'MaxNumInteractionProtein', ...
        'IndicesInteractionProtein'};

    for rv = 1:numel(requiredVars)
        if ~isfield(Sppi, requiredVars{rv})
            error('Protein file is missing required variable: %s', requiredVars{rv});
        end
    end

    A = Sppi.A;
    N = double(Sppi.N);

    ProteinLabel = Sppi.ProteinLabel;
    ProteinLabelPairs = double(Sppi.ProteinLabelPairs);

    NumInteractionProtein = double(Sppi.NumInteractionProtein);
    MaxNumInteractionProtein = double(Sppi.MaxNumInteractionProtein);
    IndicesInteractionProtein = double(Sppi.IndicesInteractionProtein);

    A_ppi = sparse(double(A));

    fprintf('PPI loaded.\n');
    fprintf('PPI size: %d x %d\n', size(A_ppi,1), size(A_ppi,2));
    fprintf('PPI edges: %d\n', nnz(triu(A_ppi,1)));
    fprintf('ProteinLabelPairs: %d x %d\n', size(ProteinLabelPairs,1), size(ProteinLabelPairs,2));

    % =========================
    % Load GO layer
    % =========================
    if ~strcmp(ExperimentType, 'TOPO')

        fprintf('\nLoading GO file:\n%s\n', GOFile);

        if ~exist(GOFile, 'file')
            error('GO file not found: %s', GOFile);
        end

        Sgo = load(GOFile);

        if ~isfield(Sgo, 'A_go')
            error('GO file does not contain variable A_go.');
        end

        A_go = sparse(double(Sgo.A_go));

        if size(A_go,1) ~= size(A_ppi,1) || size(A_go,2) ~= size(A_ppi,2)
            error('GO layer size does not match PPI adjacency size.');
        end

        A_go(A_go < 0) = 0;
        A_go(1:size(A_go,1)+1:end) = 0;

        fprintf('GO layer loaded.\n');
        fprintf('GO type: %s\n', GOType);
        fprintf('GO size: %d x %d\n', size(A_go,1), size(A_go,2));
        fprintf('GO nnz: %d\n', nnz(A_go));
        fprintf('GO max similarity: %.6f\n', full(max(A_go(:))));
        fprintf('GO diagonal nnz: %d\n', nnz(diag(A_go)));

    else

        fprintf('\nTOPO experiment selected: using empty GO layer.\n');
        A_go = sparse(size(A_ppi,1), size(A_ppi,2));

    end

    % =========================
    % Experiment flags
    % =========================
    switch ExperimentType

        case 'TOPO'
            UseGOinFitness = false;
            UseHeuristicMutation = false;
            UseRepair = false;

        case 'PPI_GO'
            UseGOinFitness = true;
            UseHeuristicMutation = false;
            UseRepair = false;

        case 'FULL'
            UseGOinFitness = true;
            UseHeuristicMutation = true;
            UseRepair = true;
    end

    fprintf('\nFlags:\n');
    fprintf('UseGOinFitness       = %d\n', UseGOinFitness);
    fprintf('UseHeuristicMutation = %d\n', UseHeuristicMutation);
    fprintf('UseRepair            = %d\n', UseRepair);

    % =========================
    % Initialize populations
    % =========================
    ResultsGroup = repmat(struct('Results', []), MaxRun, 1);
    PopulationGroup = repmat(struct('Population', []), MaxRun, 1);

    fprintf('\nInitializing populations...\n');

    for RunNumber = 1:MaxRun

        rng(InitSeeds(RunNumber), 'twister');

        fprintf('Initializing run %d / %d\n', RunNumber, MaxRun);

        PopulationGroup(RunNumber).Population = ...
            CreatePopulationOverlapping(ProteinLabelPairs, PopulationSize);
    end

    fprintf('Population initialization finished.\n');

    % =========================
    % Runtime tracking
    % =========================
    RunTimeSeconds = zeros(MaxRun,1);

    % =========================
    % Run algorithm
    % =========================
    fprintf('\nStarting optimization...\n');

    parfor RunNumber = 1:MaxRun

        rng(RunSeeds(RunNumber), 'twister');

        tRun = tic;

        fprintf('Run %d / %d started.\n', RunNumber, MaxRun);

        PopTemp = Individual2CmplxDecodingOverlapping( ...
            ProteinLabelPairs, ...
            PopulationGroup(RunNumber).Population, ...
            PopulationSize);

        PopTemp = ComputeFitnessOverlapping( ...
            A_ppi, A_go, N, ...
            IndicesInteractionProtein, ...
            NumInteractionProtein, ...
            MaxNumInteractionProtein, ...
            PopTemp, PopulationSize);

        ResultsGroup(RunNumber).Results = MOEAD_Overlapping_TNSE( ...
            A_ppi, A_go, N, ...
            IndicesInteractionProtein, ...
            NumInteractionProtein, ...
            MaxNumInteractionProtein, ...
            PopTemp, ...
            PopulationSize, ...
            Heuristic, ...
            Pm, ...
            ChromHeuristic, ...
            ProteinLabelPairs, ...
            UseGOinFitness, ...
            UseHeuristicMutation, ...
            UseRepair, ...
            Generations);

        RunTimeSeconds(RunNumber) = toc(tRun);

        fprintf('Run %d / %d finished. Time = %.2f sec\n', ...
            RunNumber, MaxRun, RunTimeSeconds(RunNumber));
    end

    TotalRuntimeSeconds = sum(RunTimeSeconds);

    fprintf('\nOptimization finished.\n');
    fprintf('Total runtime: %.2f sec\n', TotalRuntimeSeconds);
    fprintf('Mean runtime per run: %.2f sec\n', mean(RunTimeSeconds));

    % =========================
    % Save results
    % =========================
    SaveDir = fullfile('Results_TNSE', DatasetName, ExperimentType);

    if ~exist(SaveDir, 'dir')
        mkdir(SaveDir);
    end

    % File name without timestamp
    SaveName = sprintf('Results_%s_%s_Runs%d_Pop%d_Gen%d.mat', ...
        DatasetName, ExperimentType, MaxRun, PopulationSize, Generations);

    SavePath = fullfile(SaveDir, SaveName);

    save(SavePath, ...
        'ResultsGroup', ...
        'DatasetName', ...
        'ExperimentType', ...
        'NetworkNumber', ...
        'ExpChoice', ...
        'Pm', ...
        'T_go', ...
        'MaxRun', ...
        'PopulationSize', ...
        'Generations', ...
        'Heuristic', ...
        'ChromHeuristic', ...
        'UseGOinFitness', ...
        'UseHeuristicMutation', ...
        'UseRepair', ...
        'RunTimestamp', ...
        'MATLABVersion', ...
        'InitSeeds', ...
        'RunSeeds', ...
        'RunTimeSeconds', ...
        'TotalRuntimeSeconds', ...
        'GOType', ...
        'GOOntologyFile', ...
        'GOAnnotationFile', ...
        'ProteinFile', ...
        'GOFile', ...
        '-v7.3');

    fprintf('\nResults saved:\n%s\n', SavePath);
end
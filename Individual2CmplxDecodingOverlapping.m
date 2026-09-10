function [Population] = Individual2CmplxDecodingOverlapping(ProteinLabelPairs, Population, PopulationSize)
% Individual2CmplxDecodingOverlapping
% Decodes all individuals from chromosome representation to overlapping
% protein complexes.

    if nargin < 3 || isempty(PopulationSize)
        PopulationSize = numel(Population);
    end

    if PopulationSize > numel(Population)
        error('PopulationSize is larger than the actual number of individuals.');
    end

    for IndividualCounter = 1:PopulationSize

        if ~isfield(Population(IndividualCounter), 'Chromosome')
            error('Population(%d) does not contain Chromosome.', IndividualCounter);
        end

        Population(IndividualCounter).CommsNodes = ...
            ComputeCmplxDecodingOverlappingFast( ...
                Population(IndividualCounter).Chromosome, ...
                ProteinLabelPairs);
    end
end
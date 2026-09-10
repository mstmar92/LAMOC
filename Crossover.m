function [Child] = Crossover(Parent1, Parent2, N, Child)
% Crossover
% Uniform crossover for edge-based chromosome.
% N is kept for compatibility.


    Pc = 0.8;

    if ~isfield(Parent1, 'Chromosome') || ~isfield(Parent2, 'Chromosome')
        error('Both parents must contain Chromosome.');
    end

    if length(Parent1.Chromosome) ~= length(Parent2.Chromosome)
        error('Parent chromosomes have different lengths.');
    end

    numEdges = length(Parent1.Chromosome);

    mask = rand(1, numEdges) <= Pc;

    Child.Chromosome = Parent1.Chromosome;
    Child.Chromosome(~mask) = Parent2.Chromosome(~mask);
end
function [F] = ScalarFunction(Individual, Weight, IndivPoint, ObjectiveDimension)
% ScalarFunction
% Tchebycheff scalarization for minimization in MOEA/D.
%
% g(x|w,z) = max_d w_d * |f_d(x) - z_d|
%
% Obj(1): Conductance
% Obj(2): GO coherence loss

    if ~isfield(Individual, 'Obj') || numel(Individual.Obj) < ObjectiveDimension
        error('Individual does not contain a valid Obj field.');
    end

    if numel(Weight) < ObjectiveDimension
        error('Weight vector dimension is smaller than ObjectiveDimension.');
    end

    MaxFunction = -1.0e30;

    epsWeight = 1e-6;

    for d = 1:ObjectiveDimension

        if numel(IndivPoint) < d || ~isfield(IndivPoint(d), 'Ideal')
            error('IndivPoint(%d).Ideal is missing.', d);
        end

        z_d = IndivPoint(d).Ideal;
        f_d = Individual.Obj(d);

        if isnan(f_d) || isinf(f_d)
            F = Inf;
            return;
        end

        w_d = max(Weight(d), epsWeight);

        Diff = abs(f_d - z_d);
        val  = w_d * Diff;

        if val > MaxFunction
            MaxFunction = val;
        end
    end

    F = MaxFunction;
end
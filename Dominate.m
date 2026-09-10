function [NearParetoOptimalSet] = Dominate(Child, NearParetoOptimalSet, nd, ObjectiveDimension)
% Dominate
% Marks NearParetoOptimalSet(nd) as dominated if Child dominates it.
%
% Minimization dominance:
% Child dominates Old if:
%   Child.Obj(k) <= Old.Obj(k) for all objectives
%   and Child.Obj(k) < Old.Obj(k) for at least one objective.

    if ObjectiveDimension ~= 2
        error('Dominate currently supports ObjectiveDimension = 2 only.');
    end

    if ~isfield(Child, 'Obj') || numel(Child.Obj) < 2
        error('Child does not contain a valid Obj field.');
    end

    if ~isfield(NearParetoOptimalSet(nd), 'Obj') || numel(NearParetoOptimalSet(nd).Obj) < 2
        error('NearParetoOptimalSet(%d) does not contain a valid Obj field.', nd);
    end

    childObj = Child.Obj;
    oldObj   = NearParetoOptimalSet(nd).Obj;

    if any(isnan(childObj)) || any(isinf(childObj)) || ...
       any(isnan(oldObj))   || any(isinf(oldObj))

        NearParetoOptimalSet(nd).Dominated = 'F';
        return;
    end

    childDominatesOld = ...
        childObj(1) <= oldObj(1) && ...
        childObj(2) <= oldObj(2) && ...
        (childObj(1) < oldObj(1) || childObj(2) < oldObj(2));

    if childDominatesOld
        NearParetoOptimalSet(nd).Dominated = 'T';
    else
        NearParetoOptimalSet(nd).Dominated = 'F';
    end
end
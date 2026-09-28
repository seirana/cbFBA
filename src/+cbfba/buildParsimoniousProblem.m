function problem = buildParsimoniousProblem(model, complexModel, method)
%BUILDPARSIMONIOUSPROBLEM Construct cbFBA or pFBA absolute-value LP.
%
%   cbFBA minimizes sum(abs(A*v)).
%   pFBA minimizes sum(abs(v)).
%
%   Both formulations support signed/reversible reaction variables by using
%   explicit non-negative auxiliary absolute-value variables.

    model = cbfba.validateModel(model);
    method = lower(string(method));

    nReactions = numel(model.rxns);
    nMetabolites = size(model.S, 1);

    switch method
        case "cbfba"
            if nargin < 2 || isempty(complexModel) ...
                    || ~isstruct(complexModel) ...
                    || ~isfield(complexModel, "A")
                error("cbfba:MissingComplexModel", "cbFBA requires complexModel.A.");
            end
            Acomplex = sparse(complexModel.A);
            if size(Acomplex, 2) ~= nReactions
                error("cbfba:DimensionMismatch", "complexModel.A has incompatible dimensions.");
            end

            nAuxiliary = size(Acomplex, 1);
            Aineq = [ ...
                Acomplex, -speye(nAuxiliary); ...
                -Acomplex, -speye(nAuxiliary)];
            bineq = zeros(2 * nAuxiliary, 1);
            objective = [zeros(nReactions, 1); ones(nAuxiliary, 1)];

        case "pfba"
            nAuxiliary = nReactions;
            Aineq = [ ...
                speye(nReactions), -speye(nReactions); ...
                -speye(nReactions), -speye(nReactions)];
            bineq = zeros(2 * nReactions, 1);
            objective = [zeros(nReactions, 1); ones(nReactions, 1)];

        otherwise
            error("cbfba:UnknownMethod", "method must be 'cbfba' or 'pfba'.");
    end

    Aeq = [sparse(model.S), sparse(nMetabolites, nAuxiliary)];
    beq = zeros(nMetabolites, 1);
    lb = [model.lb; zeros(nAuxiliary, 1)];
    ub = [model.ub; Inf(nAuxiliary, 1)];

    problem = struct( ...
        "method", method, ...
        "objective", objective, ...
        "Aineq", Aineq, ...
        "bineq", bineq, ...
        "Aeq", Aeq, ...
        "beq", beq, ...
        "lb", lb, ...
        "ub", ub, ...
        "nReactions", nReactions, ...
        "nAuxiliary", nAuxiliary);
end

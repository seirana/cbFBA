function result = optimizeFluxRanges(model, method, toleranceFactor, options)
%OPTIMIZEFLUXRANGES Optimize cbFBA/pFBA objective and characterize alternatives.
%
%   RESULT = cbfba.optimizeFluxRanges(MODEL, METHOD, FACTOR, OPTIONS)
%   first minimizes the cbFBA or pFBA objective, then constrains the
%   objective to <= FACTOR * optimum and computes reaction-wise net flux
%   minima/maxima under that near-optimal space.
%
%   FACTOR must be >= 1. The publication examples used 1.00001.

    if nargin < 4 || isempty(options)
        options = struct();
    end
    if ~isfield(options, "complexTolerance")
        options.complexTolerance = 1e-12;
    end
    if ~isfield(options, "display")
        options.display = "none";
    end

    validateattributes(toleranceFactor, {'numeric'}, {'scalar', 'real', 'finite', '>=', 1});

    if exist("linprog", "file") ~= 2
        error("cbfba:OptimizationToolboxMissing", "Optimization Toolbox function linprog is required.");
    end

    model = cbfba.validateModel(model);
    complexModel = cbfba.buildComplexMatrices(model, options.complexTolerance);
    problem = cbfba.buildParsimoniousProblem(model, complexModel, method);

    solverOptions = optimoptions("linprog", "Display", char(string(options.display)));

    [optimumPoint, optimumValue, exitFlag, solverOutput] = linprog( ...
        problem.objective, ...
        problem.Aineq, ...
        problem.bineq, ...
        problem.Aeq, ...
        problem.beq, ...
        problem.lb, ...
        problem.ub, ...
        solverOptions);
    localRequireOptimal(exitFlag, solverOutput, "primary objective");

    if isempty(optimumPoint) || ~isfinite(optimumValue)
        error("cbfba:OptimizationFailed", "Primary optimization returned no finite optimum.");
    end

    objectiveLimit = toleranceFactor * optimumValue;
    Aineq = [problem.Aineq; problem.objective'];
    bineq = [problem.bineq; objectiveLimit];

    nReactions = problem.nReactions;
    minimum = zeros(nReactions, 1);
    maximum = zeros(nReactions, 1);

    for i = 1:nReactions
        fluxObjective = cbfba.netFluxObjective(model, i);
        objective = [fluxObjective; zeros(problem.nAuxiliary, 1)];

        [~, minimumValue, exitFlagMin, outputMin] = linprog( ...
            objective, Aineq, bineq, problem.Aeq, problem.beq, ...
            problem.lb, problem.ub, solverOptions);
        localRequireOptimal(exitFlagMin, outputMin, sprintf("minimum of %s", model.rxns{i}));

        [~, negativeMaximum, exitFlagMax, outputMax] = linprog( ...
            -objective, Aineq, bineq, problem.Aeq, problem.beq, ...
            problem.lb, problem.ub, solverOptions);
        localRequireOptimal(exitFlagMax, outputMax, sprintf("maximum of %s", model.rxns{i}));

        minimum(i) = minimumValue;
        maximum(i) = -negativeMaximum;
    end

    ranges = table( ...
        model.rxnNumber, ...
        string(model.rxns), ...
        minimum, ...
        maximum, ...
        min(abs(minimum), abs(maximum)), ...
        max(abs(minimum), abs(maximum)), ...
        'VariableNames', ...
        {'ReactionNumber', 'ReactionID', 'Minimum', 'Maximum', 'MinAbsoluteEndpoint', 'MaxAbsoluteEndpoint'});

    result = struct( ...
        "method", lower(string(method)), ...
        "optimumObjective", double(optimumValue), ...
        "objectiveLimit", double(objectiveLimit), ...
        "toleranceFactor", double(toleranceFactor), ...
        "optimumFlux", optimumPoint(1:nReactions), ...
        "ranges", ranges, ...
        "complexModel", complexModel);
end

function localRequireOptimal(exitFlag, output, context)
    if exitFlag ~= 1
        message = "";
        if isstruct(output) && isfield(output, "message")
            message = string(output.message);
        end
        error( ...
            "cbfba:OptimizationFailed", ...
            "linprog failed during %s (exitFlag=%d). %s", ...
            context, ...
            exitFlag, ...
            message);
    end
end

function model = validateModel(model)
%VALIDATEMODEL Validate and normalize a metabolic model for cbFBA.

    if ~isstruct(model)
        error("cbfba:InvalidModel", "Model must be a struct.");
    end

    required = ["S", "rxns", "mets", "lb", "ub"];
    missing = required(~isfield(model, cellstr(required)));
    if ~isempty(missing)
        error( ...
            "cbfba:MissingModelFields", ...
            "Model is missing required fields: %s", ...
            strjoin(missing, ", "));
    end

    if (~isnumeric(model.S) && ~issparse(model.S)) || ndims(model.S) ~= 2
        error("cbfba:InvalidStoichiometry", "model.S must be a two-dimensional numeric matrix.");
    end
    if any(~isfinite(nonzeros(model.S)))
        error("cbfba:InvalidStoichiometry", "model.S contains non-finite values.");
    end

    nMetabolites = size(model.S, 1);
    nReactions = size(model.S, 2);

    if numel(model.rxns) ~= nReactions || numel(model.mets) ~= nMetabolites
        error( ...
            "cbfba:DimensionMismatch", ...
            "Reaction/metabolite identifiers must match model.S dimensions.");
    end

    model.rxns = cellstr(string(model.rxns(:)));
    model.mets = cellstr(string(model.mets(:)));
    model.lb = double(model.lb(:));
    model.ub = double(model.ub(:));

    if numel(model.lb) ~= nReactions || numel(model.ub) ~= nReactions
        error("cbfba:DimensionMismatch", "lb and ub must contain one value per reaction.");
    end
    if any(~isfinite(model.lb)) || any(~isfinite(model.ub))
        error("cbfba:NonFiniteBounds", "Reaction bounds must be finite.");
    end
    if any(model.lb > model.ub)
        error("cbfba:InvalidBounds", "Every lower bound must be <= its upper bound.");
    end

    if ~isfield(model, "c") || numel(model.c) ~= nReactions
        model.c = zeros(nReactions, 1);
    else
        model.c = double(model.c(:));
    end

    if ~isfield(model, "rxnNumber") || numel(model.rxnNumber) ~= nReactions
        model.rxnNumber = (1:nReactions)';
    else
        model.rxnNumber = double(model.rxnNumber(:));
    end

    if isfield(model, "match") && ~isempty(model.match)
        match = double(model.match(:));
        if numel(match) ~= nReactions
            error("cbfba:InvalidMatchVector", "model.match must contain one value per reaction.");
        end
        if any(match < 0) || any(match > nReactions) || any(mod(match, 1) ~= 0)
            error("cbfba:InvalidMatchVector", "model.match contains invalid reaction indices.");
        end
        model.match = match;
    end
end

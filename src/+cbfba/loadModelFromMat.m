function model = loadModelFromMat(modelFile, variableName)
%LOADMODELFROMMAT Load one metabolic model struct from a MAT file safely.
%
%   MODEL = cbfba.loadModelFromMat(FILE) automatically selects the only
%   struct variable containing S, rxns, mets, lb, and ub.
%
%   MODEL = cbfba.loadModelFromMat(FILE, NAME) selects NAME explicitly.

    arguments
        modelFile (1, 1) string
        variableName (1, 1) string = ""
    end

    if ~isfile(modelFile)
        error( ...
            "cbfba:ModelFileNotFound", ...
            "Model file not found: %s", ...
            modelFile);
    end

    contents = load(modelFile);

    if strlength(variableName) > 0
        if ~isfield(contents, variableName)
            error( ...
                "cbfba:ModelVariableNotFound", ...
                "Variable %s was not found in %s.", ...
                variableName, ...
                modelFile);
        end
        candidate = contents.(variableName);
        model = cbfba.validateModel(candidate);
        return;
    end

    names = string(fieldnames(contents));
    matches = false(numel(names), 1);

    for i = 1:numel(names)
        value = contents.(names(i));
        matches(i) = isstruct(value) ...
            && all(isfield(value, {"S", "rxns", "mets", "lb", "ub"}));
    end

    candidates = names(matches);
    if isempty(candidates)
        error( ...
            "cbfba:NoModelVariable", ...
            "No metabolic-model struct was found in %s.", ...
            modelFile);
    end
    if numel(candidates) > 1
        error( ...
            "cbfba:AmbiguousModelVariable", ...
            "Multiple model-like structs were found: %s. Specify variableName.", ...
            strjoin(candidates, ", "));
    end

    model = cbfba.validateModel(contents.(candidates(1)));
end

function output = runAnalysis(model, toleranceFactor, options)
%RUNANALYSIS Validate a model and compare cbFBA with pFBA.
%
%   This is the maintained in-memory entry point. Study-specific model
%   reconstruction is intentionally kept outside the core optimization API.

    if nargin < 3 || isempty(options)
        options = struct();
    end

    output = cbfba.compareMethods( ...
        cbfba.validateModel(model), ...
        toleranceFactor, ...
        options);
end

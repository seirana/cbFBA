function paths = writeResults(outputDirectory, runName, output, metadata)
%WRITERESULTS Persist cbFBA comparison outputs after computation completes.

    arguments
        outputDirectory (1, 1) string
        runName (1, 1) string
        output struct
        metadata struct = struct()
    end

    if ~isfield(output, "comparison") ...
            || ~istable(output.comparison) ...
            || ~isfield(output, "cbFBA") ...
            || ~isfield(output, "pFBA")
        error( ...
            "cbfba:InvalidOutput", ...
            "output must be returned by cbfba.compareMethods.");
    end

    if strlength(strtrim(runName)) == 0
        error("cbfba:InvalidRunName", "runName must not be empty.");
    end

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    base = fullfile(outputDirectory, runName);
    comparisonPath = base + "_comparison.csv";
    complexPath = base + "_complexes.csv";
    metadataPath = base + "_metadata.json";
    matPath = base + "_results.mat";

    writetable(output.comparison, comparisonPath);

    complexTable = table( ...
        (1:numel(output.complexModel.complexNames))', ...
        string(output.complexModel.complexNames(:)), ...
        output.complexClassification.labels, ...
        'VariableNames', ...
        {'ComplexIndex', 'ComplexName', 'BalanceClass'});
    writetable(complexTable, complexPath);

    metadata.method = "cbFBA versus pFBA";
    metadata.toleranceFactor = output.toleranceFactor;
    metadata.cbFBAOptimumObjective = output.cbFBA.optimumObjective;
    metadata.pFBAOptimumObjective = output.pFBA.optimumObjective;
    metadata.nReactions = height(output.comparison);
    metadata.nComplexes = numel(output.complexModel.complexNames);
    metadata.generatedAtUTC = char( ...
        datetime("now", "TimeZone", "UTC", "Format", "yyyy-MM-dd'T'HH:mm:ssXXX"));

    fid = fopen(metadataPath, "w");
    if fid < 0
        error( ...
            "cbfba:FileWriteFailed", ...
            "Could not open metadata file for writing: %s", ...
            metadataPath);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, "%s", jsonencode(metadata, "PrettyPrint", true));
    clear cleanup;

    save(matPath, "output", "-v7.3");

    paths = struct( ...
        "comparison", comparisonPath, ...
        "complexes", complexPath, ...
        "metadata", metadataPath, ...
        "mat", matPath);
end

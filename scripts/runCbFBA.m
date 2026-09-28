function output = runCbFBA(modelFile, outputDirectory, toleranceFactor, variableName)
%RUNCBFBA File-driven entry point for maintained cbFBA/pFBA comparison.
%
%   OUTPUT = runCbFBA(MODELFILE, OUTPUTDIRECTORY, TOLERANCEFACTOR)
%   loads one model struct, evaluates cbFBA and pFBA using the same model and
%   near-optimal tolerance, and writes machine-readable results.
%
%   VARIABLE_NAME is optional and is only needed when MODELFILE contains
%   more than one model-like struct.

    arguments
        modelFile (1, 1) string
        outputDirectory (1, 1) string = "artifacts"
        toleranceFactor (1, 1) double {mustBeFinite, mustBeGreaterThanOrEqual(toleranceFactor, 1)} = 1.00001
        variableName (1, 1) string = ""
    end

    repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
    addpath(fullfile(repositoryRoot, "src"));

    model = cbfba.loadModelFromMat(modelFile, variableName);
    output = cbfba.runAnalysis(model, toleranceFactor, struct());

    metadata = struct( ...
        "modelFile", char(modelFile), ...
        "modelVariable", char(variableName), ...
        "matlabRelease", version("-release"), ...
        "matlabVersion", version, ...
        "gitCommit", localGitCommit(repositoryRoot), ...
        "inputSHA256", localSha256(modelFile));

    [~, runName] = fileparts(modelFile);
    output.paths = cbfba.writeResults( ...
        outputDirectory, ...
        string(runName), ...
        output, ...
        metadata);
end

function digest = localSha256(path)
    fid = fopen(path, "rb");
    if fid < 0
        error("cbfba:FileReadFailed", "Could not open %s.", path);
    end
    cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, "*uint8");
    clear cleanup;

    engine = java.security.MessageDigest.getInstance("SHA-256");
    engine.update(bytes);
    digest = lower(string(reshape(dec2hex(typecast(engine.digest(), "uint8"), 2).', 1, [])));
end

function commit = localGitCommit(repositoryRoot)
    originalDirectory = pwd;
    cleanup = onCleanup(@() cd(originalDirectory));
    cd(repositoryRoot);

    [status, textOutput] = system("git rev-parse HEAD");
    if status == 0
        commit = strtrim(string(textOutput));
    else
        commit = "";
    end
    clear cleanup;
end

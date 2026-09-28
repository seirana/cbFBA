function tests = testCbFBACore
%TESTCBFBACORE Solver-independent tests for maintained cbFBA utilities.

    tests = functiontests(localfunctions);
end

function testValidateModelAddsStableReactionNumbers(testCase)
    model = cbfba.validateModel(localToyModel());

    verifyEqual(testCase, model.rxnNumber, (1:3)');
    verifyEqual(testCase, string(model.rxns), ["R1"; "R2"; "R3"]);
end

function testValidateModelRejectsInvalidBounds(testCase)
    model = localToyModel();
    model.lb(2) = 5;
    model.ub(2) = 1;

    verifyError( ...
        testCase, ...
        @() cbfba.validateModel(model), ...
        "cbfba:InvalidBounds");
end

function testComplexMatricesReconstructStoichiometry(testCase)
    model = cbfba.validateModel(localToyModel());

    complexModel = cbfba.buildComplexMatrices(model, 1e-12);

    verifyEqual( ...
        testCase, ...
        full(complexModel.Y * complexModel.A), ...
        model.S, ...
        "AbsTol", 1e-12);
    verifySize(testCase, complexModel.A, [4, 3]);
    verifyEqual(testCase, complexModel.reconstructionError, 0, "AbsTol", 1e-12);
end

function testComplexClassificationIsDeterministic(testCase)
    complexModel = struct();
    complexModel.Y = sparse([ ...
        1, 0, 0; ...
        0, 1, 1]);
    complexModel.complexNames = ["C1"; "C2"; "C3"];

    result = cbfba.classifyComplexes(complexModel);

    verifyEqual(testCase, result.labels(1), "TB");
    verifyEqual(testCase, result.nTriviallyBalanced, 1);
    verifyEqual( ...
        testCase, ...
        result.nTriviallyBalanced ...
            + result.nNonTriviallyBalanced ...
            + result.nNonBalanced, ...
        3);
end

function testCbFBAProblemUsesComplexAbsoluteValueAuxiliaries(testCase)
    model = cbfba.validateModel(localToyModel());
    complexModel = cbfba.buildComplexMatrices(model);

    problem = cbfba.buildParsimoniousProblem( ...
        model, ...
        complexModel, ...
        "cbfba");

    verifyEqual(testCase, problem.nReactions, 3);
    verifyEqual(testCase, problem.nAuxiliary, size(complexModel.A, 1));
    verifyEqual(testCase, numel(problem.objective), 3 + problem.nAuxiliary);
    verifyEqual( ...
        testCase, ...
        problem.objective(1:3), ...
        zeros(3, 1));
    verifyEqual( ...
        testCase, ...
        problem.objective(4:end), ...
        ones(problem.nAuxiliary, 1));
end

function testPFBAProblemUsesReactionAbsoluteValueAuxiliaries(testCase)
    model = cbfba.validateModel(localToyModel());

    problem = cbfba.buildParsimoniousProblem( ...
        model, ...
        [], ...
        "pfba");

    verifyEqual(testCase, problem.nAuxiliary, 3);
    verifyEqual(testCase, size(problem.Aineq), [6, 6]);
    verifyEqual(testCase, problem.objective(1:3), zeros(3, 1));
    verifyEqual(testCase, problem.objective(4:6), ones(3, 1));
end

function testNetFluxObjectiveUsesForwardMinusReversePair(testCase)
    model = cbfba.validateModel(localToyModel());
    model.match = [3; 0; 1];

    objective = cbfba.netFluxObjective(model, 1);

    verifyEqual(testCase, objective, [1; 0; -1]);
end

function testLoadModelFromMatAutoDetectsSingleModel(testCase)
    folder = testCase.applyFixture( ...
        matlab.unittest.fixtures.TemporaryFolderFixture);
    modelPath = fullfile(folder.Folder, "toy.mat");

    toy = localToyModel(); %#ok<NASGU>
    metadata = struct("note", "not a model"); %#ok<NASGU>
    save(modelPath, "toy", "metadata");

    loaded = cbfba.loadModelFromMat(string(modelPath));

    verifyEqual(testCase, string(loaded.rxns), ["R1"; "R2"; "R3"]);
end

function model = localToyModel()
    model = struct();
    model.S = [ ...
        -1, 1, 0; ...
         0, -1, 1];
    model.rxns = {"R1"; "R2"; "R3"};
    model.mets = {"M1"; "M2"};
    model.lb = [0; -10; 0];
    model.ub = [10; 10; 10];
    model.c = zeros(3, 1);
end

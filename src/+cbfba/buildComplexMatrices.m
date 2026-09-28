function complexModel = buildComplexMatrices(model, tolerance)
%BUILDCOMPLEXMATRICES Construct the Y and incidence A matrices numerically.
%
%   For each reaction, the reactant and product stoichiometric vectors are
%   treated as chemical complexes. Identical non-zero complexes are shared
%   across reactions. Exchange-side zero complexes are omitted by default,
%   matching the historical repository behavior.

    if nargin < 2 || isempty(tolerance)
        tolerance = 1e-12;
    end
    validateattributes(tolerance, {'numeric'}, {'scalar', 'real', 'positive', 'finite'});

    model = cbfba.validateModel(model);

    S = double(model.S);
    reactants = max(-S, 0);
    products = max(S, 0);
    nReactions = size(S, 2);

    candidateComplexes = [reactants, products];
    nonzeroMask = max(abs(candidateComplexes), [], 1) > tolerance;
    nonzeroComplexes = candidateComplexes(:, nonzeroMask);

    if isempty(nonzeroComplexes)
        error("cbfba:NoComplexes", "No non-zero reaction complexes were found.");
    end

    canonical = localCanonicalize(nonzeroComplexes, tolerance);
    [~, uniqueIndices, groupIndices] = unique( ...
        canonical.', ...
        "rows", ...
        "stable");
    Y = sparse(nonzeroComplexes(:, uniqueIndices));

    candidateToGroup = zeros(1, 2 * nReactions);
    candidateToGroup(nonzeroMask) = groupIndices;

    reactantIndex = candidateToGroup(1:nReactions);
    productIndex = candidateToGroup(nReactions + 1:end);

    nComplexes = size(Y, 2);
    A = sparse(nComplexes, nReactions);

    for reactionIndex = 1:nReactions
        if reactantIndex(reactionIndex) > 0
            A(reactantIndex(reactionIndex), reactionIndex) = -1;
        end
        if productIndex(reactionIndex) > 0
            A(productIndex(reactionIndex), reactionIndex) = 1;
        end
    end

    complexNames = strings(nComplexes, 1);
    for complexIndex = 1:nComplexes
        complexNames(complexIndex) = localComplexName( ...
            Y(:, complexIndex), ...
            model.mets);
    end

    reconstructionError = norm(model.S - Y * A, "fro");
    scale = max(1, norm(model.S, "fro"));
    if reconstructionError > tolerance * scale * 10
        error( ...
            "cbfba:ComplexDecompositionFailed", ...
            "Y*A does not reconstruct S within tolerance (error=%g).", ...
            reconstructionError);
    end

    complexModel = struct( ...
        "Y", Y, ...
        "A", A, ...
        "complexNames", complexNames, ...
        "reactantComplexIndex", reactantIndex(:), ...
        "productComplexIndex", productIndex(:), ...
        "reconstructionError", reconstructionError, ...
        "tolerance", tolerance);
end

function canonical = localCanonicalize(values, tolerance)
    scale = max(1, max(abs(values), [], "all"));
    digits = max(0, ceil(-log10(tolerance / scale)));
    digits = min(digits, 14);
    canonical = round(values, digits);
end

function name = localComplexName(column, metaboliteIDs)
    indices = find(abs(column) > 0);
    pieces = strings(numel(indices), 1);

    for i = 1:numel(indices)
        index = indices(i);
        pieces(i) = string(full(column(index))) + "*" + string(metaboliteIDs{index});
    end

    name = strjoin(pieces, "+");
end

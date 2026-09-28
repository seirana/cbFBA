function result = classifyComplexes(complexModel)
%CLASSIFYCOMPLEXES Classify complexes with deterministic leaf elimination.
%
%   RESULT = cbfba.classifyComplexes(COMPLEXMODEL) classifies complexes as:
%     "TB"  - trivially balanced by a metabolite unique to that complex;
%     "NTB" - becomes uniquely identifiable after previously balanced
%             complexes are removed;
%     "NB"  - remains non-balanced under this structural criterion.
%
%   This routine formalizes the iterative intent of the historical TBC/NTBC
%   scripts without mutating the matrix while iterating over rows.

    if ~isstruct(complexModel) ...
            || ~isfield(complexModel, "Y") ...
            || ~isfield(complexModel, "complexNames")
        error( ...
            "cbfba:InvalidComplexModel", ...
            "complexModel must contain Y and complexNames.");
    end

    Y = sparse(complexModel.Y);
    nComplexes = size(Y, 2);

    if numel(complexModel.complexNames) ~= nComplexes
        error( ...
            "cbfba:DimensionMismatch", ...
            "complexNames must contain one name per column of Y.");
    end

    labels = repmat("NB", nComplexes, 1);
    active = true(nComplexes, 1);

    initial = localLeafComplexes(Y, active);
    labels(initial) = "TB";
    active(initial) = false;

    changed = true;
    while changed
        next = localLeafComplexes(Y, active);
        next = next(active(next));

        if isempty(next)
            changed = false;
        else
            labels(next) = "NTB";
            active(next) = false;
        end
    end

    result = struct( ...
        "labels", labels, ...
        "complexNames", string(complexModel.complexNames(:)), ...
        "nTriviallyBalanced", nnz(labels == "TB"), ...
        "nNonTriviallyBalanced", nnz(labels == "NTB"), ...
        "nNonBalanced", nnz(labels == "NB"));
end

function indices = localLeafComplexes(Y, active)
    if ~any(active)
        indices = zeros(0, 1);
        return;
    end

    incidence = abs(Y(:, active)) > 0;
    occurrenceCount = full(sum(incidence, 2));
    leafMetabolites = find(occurrenceCount == 1);

    if isempty(leafMetabolites)
        indices = zeros(0, 1);
        return;
    end

    activeIndices = find(active);
    localIndices = zeros(numel(leafMetabolites), 1);

    for i = 1:numel(leafMetabolites)
        localIndices(i) = find( ...
            incidence(leafMetabolites(i), :), ...
            1, ...
            "first");
    end

    indices = unique(activeIndices(localIndices), "stable");
end

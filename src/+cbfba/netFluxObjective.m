function objective = netFluxObjective(model, reactionIndex)
%NETFLUXOBJECTIVE Build the historical net-flux objective for one reaction.
%
%   For an irreversible split pair, when model.match(i) > i the forward
%   reaction is combined as v_i - v_match(i). Otherwise the objective is v_i.

    model = cbfba.validateModel(model);
    validateattributes(reactionIndex, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', numel(model.rxns)});

    objective = zeros(numel(model.rxns), 1);
    objective(reactionIndex) = 1;

    if isfield(model, "match") && ~isempty(model.match)
        partner = model.match(reactionIndex);
        if partner > reactionIndex
            objective(partner) = -1;
        end
    end
end

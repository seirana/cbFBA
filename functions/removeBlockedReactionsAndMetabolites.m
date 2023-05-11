function model = removeBlockedReactionsAndMetabolites(model)
%remove blocked reactions and metabolites, zero columns and metabolites
%
% USAGE:
%     model = removeBlockedReactionsAndMetabolites(model)
%
% INPUT:
%     model: the metabolic network with fields:
%         * .S -       the associated sparse stoichiometric matrix
%         * .rxns -    the cell array of reactions abbreviations
%         * .mets -    the cell array of metabolite abbreviations
%         * .lb -      the doulbe array of reaction flux lower bound
%         * .ub -      the doulbe array of reaction flux upper bound
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .S -       the associated sparse stoichiometric matrix
%         * .rxns -    the cell array of reactions abbreviations
%         * .mets -    the cell array of metabolite abbreviations
%         * .lb -      the doulbe array of reaction flux lower bound
%         * .ub -      the doulbe array of reaction flux upper bound

%remove blocked reactions from model
sz = size(model.rxns,1);
for i = sz:-1:1
    sm = 0;
    for j = 1:size(model.mets,1)
        sm = sm + abs(model.S(j,i));
    end
    if sm == 0 || (model.lb(i) == 0 && model.ub(i) == 0)
        model.S(:,i) = 0;
        model.lb(i) = 0;
        model.ub(i) = 0;
    end
end

%remove dead-end metabolites from model
sz = size(model.mets,1);
for i = sz:-1:1
    sm = 0;
    for j = 1:size(model.rxns,1)
        sm = sm + abs(model.S(i,j));
    end
    if sm == 0
        model.S(i,:) = 0;
    end    
end

end
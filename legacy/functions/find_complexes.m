function model = find_complexes(model)
% remove blocked reactions and metabolites, zero columns and metabolites
% finds complexes of the model
% builds Y matrix
% builds A matrix
% finds trivially balnced complexes
% finds non-trivially balanced complexes
% finds number of balnced complexes
%
% USAGE:
%     model = find_complexes(model)
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
%         * .cmpxName -  the cell array of complxes names
%         * .Y -  the associated sparse Y matrix
%         * .A -  the associated sparse A matrix
%         * .complex_type - the type of complexes, balanced, trivially balanced, or non-balanced 
%         * .ingoingRxnName - the cell array contains incoming reacion names per complex
%         * .outgoingRxnName - the cell array contains outgoing reacion names per complex
%         * in - numner of incoming reacions per complex
%         * .out - numner of outgoing reacions per complex
%

model = removeBlockedReactionsAndMetabolites(model); % remove blocked reactions and dead-end metabolited
model = complexes(model); % find complexes
model = Y_matrix(model); % build Y matrix
model = A_matrix(model); % build A matrix
model = TBC(model); % find trivially balanced complexes
model = NTBC(model); %find non-trivially balanced complexes
model = in_out_rxns(model); % in/out reactions per complexes

end
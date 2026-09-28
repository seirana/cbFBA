function model = in_out_rxns(model)
% find incoming/outgoing reactions to trivially balanced complexes
%
% USAGE:
%     model = in_out_rxns(model)
%
% INPUT:
%     model: the metabolic network with fields:
%         * .A - the associated sparse A matrix
%         * .rxns -    the cell array of ractions abbreviations
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .ingoingRxnName - the cell array contains incoming reacion names per complex
%         * .outgoingRxnName - the cell array contains outgoing reacion names per complex
%         * in - numner of incoming reacions per complex
%         * .out - numner of outgoing reacions per complex

model.in = zeros(size(model.A,1),1);
model.out = zeros(size(model.A,1),1);

for i = 1:size(model.A,1)
    model.incomingRxnName(i,1:size(string(model.rxns(model.A(i,:) > 0,1)),1)) = string(model.rxns(model.A(i,:) > 0,1));
    model.in(i,1) = size(string(model.rxns(model.A(i,:) > 0,1)),1);
       
    model.outgoingRxnName(i,1:size(string(model.rxns(model.A(i,:) < 0,1)),1)) = string(model.rxns(model.A(i,:) < 0,1));
    model.out(i,1) = size(string(model.rxns(model.A(i,:) < 0,1)),1);
end
end
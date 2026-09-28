function model = Y_matrix(model)
% builds Y matrix
%
% USAGE:
%     model = Y_matrix(model)
%
% INPUT:
%     model: the metabolic network with fields:
%         * .S -       the associated sparse stoichiometric matrix
%         * .mets -    the cell array of metabolite abbreviations
%         * .cmpxName -  the cell array of complxes names
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .Y -  the associated sparse Y matrix


model.Y = zeros(size(model.S,1),size(model.cmpxName,1));

for i = 1:size(model.cmpxName,1)
    sep_cmpx = split(string(model.cmpxName(i,1)),'+');
    for j = 1:size(sep_cmpx,1)
        a = split(string(sep_cmpx(j,1)),'*');
        model.Y(model.mets == a(2,1),i) = str2double(a(1,1));
    end
end

end
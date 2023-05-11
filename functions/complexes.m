function model = complexes(model)
% finds the complexes of the model
%
% USAGE:
%     model = complexes(model)
%
% INPUT:
%     model: the metabolic network with fields:
%         * .mets -    the cell array of metabolite abbreviations
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .cmpxName -  the cell array of complxes names

str_S = string(abs(model.S));

model.cmpxName = cell(1,1);%complexes
c = 0;
for r = 1:size(model.S,2)
    c = c+1;
    model.cmpxName(c,1) = {''};
    for m = 1:size(model.S,1)
        if model.S(m,r) < 0
            if strcmp(model.cmpxName(c,1),'') == 1
                model.cmpxName(c,1) = {strcat(str_S(m,r),'*',string(model.mets(m,1)))};
            else
                model.cmpxName(c,1) = {strcat(string(model.cmpxName(c,1)),'+', str_S(m,r),'*',string(model.mets(m,1)))};
            end
        end
    end
    if strcmp(string(model.cmpxName(c,1)), '') == 1
        c = c-1;
    end
    
    c = c+1;
    model.cmpxName(c,1) = {''};
    for m = 1:size(model.S,1)
        if model.S(m,r) > 0
            if strcmp(model.cmpxName(c,1),'') == 1
                model.cmpxName(c,1) = {strcat(str_S(m,r),'*',string(model.mets(m,1)))};
            else
                model.cmpxName(c,1) = {strcat(string(model.cmpxName(c,1)),'+', str_S(m,r),'*',string(model.mets(m,1)))};
            end
        end
    end
    if strcmp(string(model.cmpxName(c,1)), '') == 1
        c = c-1;
        
        model.cmpxName(end,:) = [];
    end
end
model.cmpxName = unique(string(model.cmpxName),'stable');

end
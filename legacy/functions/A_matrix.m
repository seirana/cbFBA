function model = A_matrix(model)
% builds A matrix
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
%         * .A -  the associated sparse A matrix


str_S = string(abs(model.S));

%build matrix A
model.A = zeros(size(model.cmpxName,1),size(model.S,2));

for r = 1:size(model.S,2)
    tmp(1,1) = {''};
    
    for m = 1:size(model.S,1)
        if model.S(m,r) < 0
            if strcmp(tmp(1,1),'') == 1
                tmp(1,1) = {strcat(str_S(m,r),'*',string(model.mets(m,1)))};                           
            else
                tmp(1,1) = {strcat(string(tmp(1,1)),'+', str_S(m,r),'*',string(model.mets(m,1)))};
            end
        end
    end
    model.A(string(model.cmpxName) == string(tmp(1,1)),r) = -1;
    
    tmp(1,1) = {''};
    for m = 1:size(model.S,1)
        if model.S(m,r) > 0
            if strcmp(tmp(1,1),'') == 1
                tmp(1,1) = {strcat(str_S(m,r),'*',string(model.mets(m,1)))};
            else
                tmp(1,1) = {strcat(string(tmp(1,1)),'+', str_S(m,r),'*',string(model.mets(m,1)))};
            end
        end
    end
    model.A(string(model.cmpxName) == string(tmp(1,1)),r) = 1;
end

end
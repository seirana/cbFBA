function model = TBC(model)
% finds trivially balnced complexes
%
% USAGE:
%     model = TBC(model)
%
% INPUT:
%     model: the metabolic network with fields:
%         * .cmpxName - the cell array of complxes names
%         * .Y - the associated sparse Y matrix
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .complex_type - the type of complexes, balanced, trivially balanced, or non-balanced 

model.complex_type(1:size(model.cmpxName,1),1) = "non-balanced";

for i = 1:size(model.Y,1)
    tc = find(model.Y(i,:) ~= 0);
    if size(tc,2) == 1
        model.complex_type(tc,1) = "TB";  
    end
end

end
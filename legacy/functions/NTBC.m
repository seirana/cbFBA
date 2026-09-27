function model = NTBC(model)
% finds non-trivially balanced complexes
%
% USAGE:
%     model = TBC(model)
%
% INPUT:
%     model: the metabolic network with fields:
%         * .Y - the associated sparse Y matrix
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .complex_type - the type of complexes, trivially balanced, non-trivially balanced, or non-balanced

Y_nontrv = model.Y;
br = 0;
while br == 0
    tr = model.complex_type ~= "non-balanced";
    Y_nontrv(:,tr) = 0;
    for i = 1:size(Y_nontrv,1)
        tc = find(Y_nontrv(i,:) ~= 0);
        if size(tc,2) == 1
            model.complex_type(tc,1) = "NTB";
        else
            br = 1;
            break;
        end
    end
end

end
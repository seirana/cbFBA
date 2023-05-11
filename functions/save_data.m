function save_data(sub_model, netfluxrange, i)
% saves the results
%
% USAGE:
%     save_data(sub_model, netfluxrange, i)
%
% INPUT:
%    sub_model: the methabolic network
%    metfluxrange: calculated net flux ranges with cbFBA and pFBA methods
%    i: counter for the sub model
%

%title of saved net flux ranges
title = ["rxns","cbFBA","",  "pFBA","",;...
    "",    "lb",   "ub","lb",  "ub"];
netfluxrange = [string(sub_model.rxns),netfluxrange];
netfluxrange = [title;netfluxrange];

%save net flux ranges
writematrix(netfluxrange,strcat('netfluxrange_',string(i),'.xlsx'));

% save the sub model
save(strcat('submodel_',string(i),'.mat'),'sub_model');
end

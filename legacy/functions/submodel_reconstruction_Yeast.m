function  sub_model = submodel_reconstruction_Yeast(c)
%reconstruction sub models baes on Blank et al study

[~,~,CS] = xlsread('Yeast_CarbonSources.xlsx');
CS = string(CS(2:end,:));

[~,~,conds] = xlsread('Yeast.xlsx','A');
conds = string(conds);

[~,~,comp] = xlsread('Yeast.xlsx','REACTIONS');
comp = string(comp(3:end,14:end));

load('Yeast.mat');

for i = 1:size(model.ub,1)
    if model.ub(i,1) > 100
        model.ub(i,1) = 100;
    end
    if model.lb(i,1) < -100
        model.lb(i,1) = -100;
    end    
end

rnd = 0.000001;

% deactive all carbon sources / glucose
for t = 1:size(CS,1)
    model.lb(str2double(CS(t,1))) = 0;
end
% active glocuse
model.lb(str2double(CS(3,1))) = -100; %'r_1714' = 'D-glucose exchange'

f = zeros(7,1);
for i = 1:7
    f(i,1) = find(model.rxnNames == conds(1,2*i));
end

irr_model = convertToIrreversible(model);%model to irr.model


lu = str2double(conds(c,2:end-1));
for i = 1:7
    if lu(1,2*i) < rnd
        lu(1,2*i) = rnd;
    end
    if lu(1,2*i-1) > 0
        
        irr_model.lb(f(i,1),1) = lu(1,2*i-1)-lu(1,2*i);
        
        if lu(1,2*i-1)-lu(1,2*i) < 0
            irr_model.lb(f(i,1),1) = 0;
        end
        irr_model.ub(f(i,1),1) = lu(1,2*i-1)+lu(1,2*i);
        
    else

        irr_model.lb(irr_model.match(f(i,1),1),1) = -lu(1,2*i-1)-lu(1,2*i);
        
        if -lu(1,2*i-1)-lu(1,2*i) < 0
            irr_model.lb(irr_model.match(f(i,1),1),1) = 0;
        end
        irr_model.ub(irr_model.match(f(i,1),1),1) = -lu(1,2*i-1)+lu(1,2*i);
    end    
end
end
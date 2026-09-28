function  sub_model = submodel_reconstruction_Ecoli(i)
%reconstruction of submodels based on McCloskey et al 2018 studies

load('Ecoli.mat'); % load model
load('Ecoli_strains'); % load strains
[~,~,D_CS] = xlsread('Ecoli_DavidiCarbonSources.xlsx'); % load carbon sources
[~,~,lab] = xlsread('Ecoli_experimental_data.xlsx'); % load experimental data
lab_data = lab(4:end,2:end);
lab_data = str2double(string(lab_data));
relaxation = 0;
eps = 0.0001;
[irr_model, KO_lb, KO_ub] = parse_KO_data(iJO1366, 'core', KOdata_file, lab_data, D_CS, relaxation,eps); % build up submodels

%upate the flux ranges of the strains
sub_model = irr_model;
for j = 1:3229
    if KO_ub(j,i) > 0
        sub_model.lb(j,1) = KO_lb(j,i);
        sub_model.ub(j,1) = KO_ub(j,i);
    end
end

end
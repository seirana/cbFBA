function [irr_model, KO_lb,KO_ub] = parse_KO_data(model, biomass_choose, KOdata_file, D_CS, relaxation, eps)
% parse_KO_data generates the strains "PLZ check the paper for references"
%
% USAGE:
%     [irr_model, KO_lb,KO_ub] = parse_KO_data(model, biomass_choose, KOdata_file, D_CS, relaxation, eps)
%
% INPUT:
%     model: the metabolic network with fields:
%         * .S -       the associated sparse stoichiometric matrix
%         * .rev -     the 0-1 indicator vector of the reversible reactions
%         * .genes -   the cell array with the name of the genes
%         * .grRules - the cell representation of the GPR rules defined in a readable format.
%         * .mets -    the cell array of metabolite abbreviations
%         * .lb -      the doulbe array of reaction flux lower bound
%         * .ub -      the doulbe array of reaction flux upper bound
%
%     biomass_choose: select the active biomass reaction, 'core' for core biomass or 'wt' for wild type 
%     KOdata_file: read the information to build the strains
%     D_CS: carbon sources
%     relaxation: a vaue to active/ dective some reactions
%     eps: error rate for the fluxes
%
% OUTPUT:
%     irr_model: the metabolic network with fields:
%         * .S -         the associated sparse stoichiometric matrix
%         * .rev -       the 0-1 indicator vector of the reversible reactions
%         * .genes -     the cell array with the name of the genes
%         * .grRules -   the cell representation of the GPR rules defined in a readable format.
%         * .mets -      the cell array of metabolite abbreviations
%         * .lb -        the doulbe array of reaction flux lower bound
%         * .ub -        the doulbe array of reaction flux upper bound
%         * .rxnNumber - the number of reactions in the original model
%             
%     KO_lb: lower bound for some selected reaction per strain  
%     KO_ub: upper bound for some selected reaction per strain 

%remove the s0001 grRules for the model. There are 3 types
%'s0001', 's0001 or ...' and 'bxxx or s0001'

excep=find(contains(model.grRules,'s0001'));
copy=model.grRules(excep);
copy=erase(copy,'s0001 or ');
copy=erase(copy,' or s0001');
copy=erase(copy,'s0001');
model.grRules(excep)=copy;

%remove one biomass reaction from model
if strcmp(biomass_choose,'wt')
    bio = find(contains(model.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M'));
    model.lb(bio) = 0;
    model.ub(bio) = 1000;
    
    f = find(contains(model.rxns,'BIOMASS_Ec_iJO1366_core_53p95M'));
    %model=removeRxns(model,model.rxns(f));
    model.lb(f) = 0;
    model.ub(f) = 0;
    model.S(:,f) = 0;
end

if strcmp(biomass_choose,'core')
    bio = find(contains(model.rxns,'BIOMASS_Ec_iJO1366_core_53p95M'));
    model.lb(bio) = 0;
    model.ub(bio) = 1000;
    
    f = find(contains(model.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M'));
    %model=removeRxns(model,model.rxns(14));
    model.lb(f) = 0;
    model.ub(f) = 0;
    model.S(:,f) = 0;
end

[~,cs_ind] = ismember(D_CS, model.rxns);%active DAVIDI exchange carbon sources
model.lb(cs_ind) = -1000;

irr_model = convertToIrreversible(model);%model to irr. model

%add rxn numbers and rxn type to the model
irr_model.rxnNumber = zeros(size(irr_model.rxns,1),1);
for j = 1:size(irr_model.rxns,1)
    irr_model.rxnNumber(j,1) = j;
end

D_CS = append(D_CS,'_b');
[~,cs_ind] = ismember(D_CS, irr_model.rxns);%deactive all DAVIDI exchange carbon sources
irr_model.ub(cs_ind) = 0;

%create lb and ub matrices
n_r = size(irr_model.S,2);
KO_lb = zeros(n_r,18);
KO_ub = zeros(n_r,18);
for i = 1:18
    KO_lb(:,i) = irr_model.lb;
    KO_ub(:,i) = irr_model.ub;
    
    KO_genes = KOdata_file.strain_mutations{i};
    [~,gene_ind] = ismember(KO_genes,irr_model.proteins);
    gene_ind = gene_ind(gene_ind~= 0);
    
    reac_ind = find(contains(irr_model.grRules,irr_model.genes(gene_ind)));
    
    %for every reaction containing at least one of theses genes
    for j = 1:length(reac_ind)
        check = split(irr_model.grRules(reac_ind(j)),' or ');
        check2 = split(irr_model.grRules(reac_ind(j)),' and ');
        if length(check)>1 && length(check2) == 1%isoenzyme
            
            [~,ind2] = ismember(check,irr_model.genes(gene_ind));%check if the genes belong to the KO set
            if sum(ind2~= 0) == length(check)%all possible isoenzymes are KO
                KO_ub(reac_ind(j),i) = relaxation;
                
            end
            
        elseif length(check2) > 1 && length(check) == 1%only having 'and'
            KO_ub(reac_ind(j),i) = relaxation;
            
            
        elseif length(check) == 1 && length(check2) == 1%single genes
            
            KO_ub(reac_ind(j),i) = relaxation;
            
        elseif length(check) > 1 && length(check2)>1%complex (&)or(&)
            count = 0;
            for m = 1:length(check)
                check3 = split(check(m),' and ');
                check3 = erase(check3,'(');
                check3 = erase(check3,')');
                
                [~,ind] = ismember(check3,irr_model.genes(gene_ind));
                if sum(ind ~= 0) > 0
                    count = count+1;
                end
            end
            if count == length(check)
                KO_ub(reac_ind(j),i) = relaxation;
                
            end
        end
    end
    
    %assign value for biomass (lb and ub)
    if KOdata_file.gro_rate(i,2) ~= 0
        KO_lb(bio,i) = KOdata_file.gro_rate(i,1)-eps;
        KO_ub(bio,i) = KOdata_file.gro_rate(i,2)+eps;
    end
    
    %assign value for Glucose uptake rxn
    f = find(contains(irr_model.rxns,'EX_glc__D_e_b'));
    if KOdata_file.glc_up(i,2) ~= 0
        KO_lb(f,i) = KOdata_file.glc_up(i,1)-eps;
        KO_ub(f,i) = KOdata_file.glc_up(i,2)+eps;
    end
    
    %assign value for Acetate secration rxn
    f = find(contains(irr_model.rxns,'EX_ac_e_f'));
    if KOdata_file.ace_se(i,2) ~= 0
        KO_lb(f,i) = KOdata_file.ace_se(i,1)-eps;
        KO_ub(f,i) = KOdata_file.ace_se(i,2)+eps;
    end
    
    %assign value for Succinate sexration rxn
    f = find(contains(irr_model.rxns,'EX_succ_e_f'));
    if KOdata_file.suc_se(i,1) ~= 0
        KO_lb(f,i) = KOdata_file.suc_se(i,1)-eps;
        KO_ub(f,i) = KOdata_file.suc_se(i,2)+eps;
    end
end
end
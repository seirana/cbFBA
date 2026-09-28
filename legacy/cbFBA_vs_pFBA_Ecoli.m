clc;clear;

err = 1.00001;
for i = 1:17        
    sub_model = reconstruction_submodels_Ecoli(i);
    sub_model = find_complexes(sub_model);    
    
    a = cbFBA(sub_model,err);
    netfluxrange(:,1:2) = a;
    b = pFBA(sub_model,err);
    netfluxrange(:,3:4) = b;
    
    save_data(sub_model, netfluxrange);
end
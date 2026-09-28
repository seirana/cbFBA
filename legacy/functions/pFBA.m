function lb_ub = pFBA(IM,err)
% calculates the flux ranges with pFBA method
%
% USAGE:
%     lb_ub = pFBA(IM, err)
%
% INPUT:
%    IM: the methabolic network with fields:
%         * .S -   the associated sparse stoichiometric matrix
%         * .rxns -    the cell array of ractions abbreviations
%         * .cmpxName -  the cell array of complxes names
%         * .A -  the associated sparse A matrix
%         * .lb -      the doulbe array of reaction flux lower bound
%         * .ub -      the doulbe array of reaction flux upper bound
%         * .match -      the double array shows if reaction is
%         reversible (= reverse reaction number) or irrevesible( =0)
%
% OUTPUT:
%     lb_ub = net flux range for reactions after applying pFBA  

options = optimset('linprog');
options.Display = 'off';

c = size(IM.cmpxName,1);
r = size(IM.S,2);
m = size(IM.S,1);

bk = endsWith(IM.rxns,'_b');
for b = 1:size(bk,1)
    if bk(b) == 1
        break;
    end
end
b = b-1;

fk = [ones(r,1);zeros(c,1)];%factors for objective function

Ain1 = [IM.A,-1*eye(c)];%c <= t
bin1 = zeros(c,1);
Ain2 = [-1*IM.A,-1*eye(c)];%-t <= c
bin2 = zeros(c,1);

Ain = [Ain1;Ain2];
bin = [bin1;bin2];

Aeq = [IM.S,zeros(m,c)];%Nv = YAv = 0
beq = zeros(size(Aeq,1),1);

lb = [IM.lb;zeros(c,1)];%lb
ub = [IM.ub;Inf(c,1)];%ub

%find min sum flux
minFlux = linprog(fk,Ain,bin,Aeq,beq,lb,ub,options);
minFlux = round(minFlux,6);

lb_ub = zeros(r,2);

if size(minFlux,1) > 0
    minFluxSum = sum(minFlux(1:r,1));%fk = sum of rxn fluxes
    
    %add a new constraint to the model
    Ain = [Ain;fk'];
    bin = [bin;err*minFluxSum];
    
    %FVA
    for i = 1:r
        IM.c = zeros(r,1);
        IM.c(i,1)=1;
        if IM.match(i,1) > i
            IM.c(IM.match(i,1))=-1;
        end
        
        fk = [IM.c;zeros(c,1)];%factors for objective function
        [~,Sol.f,Sol.stat,~] = linprog(-1*fk,Ain,bin,Aeq,beq,lb,ub,options);
        if Sol.stat~=1
            l=0;
        else
            l=Sol.f*-1;
        end
        [~,Sol.f,Sol.stat,~] = linprog(fk,Ain,bin,Aeq,beq,lb,ub,options);
        if Sol.stat~=1
            u=0;
        else
            u=Sol.f;
        end
        if abs(u) >= abs(l)
            lb_ub(i,1) = abs(l);
            lb_ub(i,2) = abs(u);
        else
            lb_ub(i,1) = abs(u);
            lb_ub(i,2) = abs(l);
        end
    end
end
end
function random_flux_sampling(IM,n)
%make a random set of points from the feasible space
%
% USAGE:
%     random_flux_sampling(model,n)
%
% INPUT:
%    n: number of random samples
%    IM: the methabolic network with fields:
%         * .S -   the associated sparse stoichiometric matrix
%         * .lb -      the doulbe array of reaction flux lower bound
%         * .ub -      the doulbe array of reaction flux upper bound
%


options = optimset('quadprog','algorithm','active-set');
options.Display = 'off';

r = size(IM.S,2);
m = size(IM.S,1);

H = 2*eye(r);

A = [];
b = [];

Aeq = IM.S;%Nv = YAv = 0
beq = zeros(size(Aeq,1),1);

lb = IM.lb;%lb
ub = IM.ub;%ub

%find random samples
for i = 1:n
    
    rnd = rand(r,1);    
    rnd = (ub - lb) * rnd+ lb;
    
    f = -2*rnd;
    x = quadruple(H,f, A, b, Aeq, beq, lb, ub, [], options);
    sample_set(:,end*1) = round(x,6);
end

% save the sample set
save('sample_set.mat','sample_set');
end
%% ========================================================================
%  score function (Multivariate log-normal) 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [output]=scfcn_lgnrm(X, mu, par_mu, cova)
nd = 10;
% mu: paramter of lognormal 
% cova: paramter of lognormal 
DM = zeros(nd, nd);
for i=1:nd 
    DM(i,i) = (mu(i)^2 + 2*0.05^2)/(mu(i)^3 + mu(i)*0.05^2);
end 

% X1-X10
dist = log(X(1:nd)) - par_mu; 

output = DM*(cova\dist');

end 
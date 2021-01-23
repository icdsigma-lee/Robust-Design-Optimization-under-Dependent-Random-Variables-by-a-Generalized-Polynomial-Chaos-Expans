%% ========================================================================
% Multivariate normal PDF  
% Input: mean vector and covariance matrix of X 
% output: multivariate normal PDF value 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================

function [output]=nGaussDen(X, mu, cov)

N = length(X);
mahDist = (X-mu)'/cov*(X-mu);
detCov = det(cov);
num = exp((-1/2)*mahDist); 
denom = sqrt((2*pi)^N*detCov);
output = num/denom; 

end 
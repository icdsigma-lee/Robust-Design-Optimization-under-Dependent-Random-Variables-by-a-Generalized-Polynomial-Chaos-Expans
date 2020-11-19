%% Multivariate Gauss distribution density function 
% The code gets input as X, mean row-vector and covariance matrix 
% outputs multivariate Gaussian PDF value 
% The code is written by Dongjin Lee 6/19/2018
function [output]=nGaussDen(X, mu, cov)

N = length(X);
   
mahDist = (X-mu)'/cov*(X-mu);
detCov = det(cov);
num = exp((-1/2)*mahDist); 
denom = sqrt((2*pi)^N*detCov);
output = num/denom; 

end 
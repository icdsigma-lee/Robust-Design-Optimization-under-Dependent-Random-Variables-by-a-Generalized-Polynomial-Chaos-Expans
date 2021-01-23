%% Multivariate Gauss distribution density function 
% The code gets input as X, mean row-vector and covariance matrix 
% outputs multivariate Gaussian PDF value 
% The code is written by Dongjin Lee 6/19/2018
function [output]=nGaussDen(X, mu, cov)

N = length(X);
% Example of mean vector and covariance matrix 
%mu1 = 1; mu2 = 1;
%sig1 = 1; sig2 = 1;
%cov12 = 0.12;
%mu = [mu1; mu2];
%cov = [sig1^2 cov12*sig1*sig2; cov12*sig1*sig2 sig2^2];

% sqrt(mahdist) indicates Mahalanobis distance 
% if N == 1 
%     output = (1/(sqrt(2*pi)*cov))*exp(-(X-mu)^2/(2*cov^2));
% else 
   
mahDist = (X-mu)'/cov*(X-mu);
detCov = det(cov);
num = exp((-1/2)*mahDist); 
denom = sqrt((2*pi)^N*detCov);
output = num/denom; 

end 
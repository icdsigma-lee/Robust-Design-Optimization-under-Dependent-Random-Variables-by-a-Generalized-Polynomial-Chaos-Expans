%% ========================================================================
%  score function (Multivariate Gaussian) 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [output]=scoreGauss(X, mu, cov)
%Multivariate Gaussian score function  

dist = mu.*X-mu;
output = cov\dist';

end 
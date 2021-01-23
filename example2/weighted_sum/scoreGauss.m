%% ========================================================================
%  score function (Multivariate Gaussian) 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [output]=scoreGauss(X, mu, cov)
%Multivariate Gaussian score function  
% shift transformation remain dist x the same 
dist = X;
output = cov\dist';

end 
%% ========================================================================
%  y1 function (example1-2)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [output]=responY2(X,mu)
global cntY2
cntY2 = cntY2 + 1;
%Multivariate Gaussian score function  
x1 = mu(1)*X(1);
x2 = mu(2)*X(2);
output = x1 + x2 - 6.45;
end 
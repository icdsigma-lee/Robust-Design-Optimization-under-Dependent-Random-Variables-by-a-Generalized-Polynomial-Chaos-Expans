%% ========================================================================
%  y0 function (example1)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [output]=responY1(X,mu)
global cntY1
cntY1 = cntY1 + 1;
%Multivariate Gaussian score function  
x1 = mu(1)*X(1);
x2 = mu(2)*X(2);
output = (x1-4)^3 + (x1-3)^4 + (x2-5)^2+10;
end 
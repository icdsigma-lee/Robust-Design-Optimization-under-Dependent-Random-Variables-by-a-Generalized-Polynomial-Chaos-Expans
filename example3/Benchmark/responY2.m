%% ========================================================================
%  y2 function (example2)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [output]=responY2(X,mu)
%Multivariate Gaussian score function  
global cntY2;
cntY2 = cntY2 + 1;
x1 = X(1)*mu(1)*1E-4;
x2 = X(2)*mu(2)*1E-4;
x3 = X(3)*mu(3);
x4 = X(4)*mu(4);
x5 = X(5)*mu(5);
x6 = X(6)*mu(6);
x7 = X(7)*mu(7);

output = (10*x7*sqrt(1 + x4^2))*(8*x3 - 1)/(sqrt(65)*x2*x6*(x3 + x4)) - 1;

end 
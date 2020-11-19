%% ========================================================================
%  y0 function (example2)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [output]=responY0(X,mu)
global cntY0
cntY0 = cntY0 + 1;

x1 = X(1)*mu(1)*1E-4;
x2 = X(2)*mu(2)*1E-4;
x3 = X(3)*mu(3);
x4 = X(4)*mu(4);
x5 = X(5)*mu(5);

output = x5*(x1*sqrt(1+x3^2) + x2*sqrt(1+x4^2));

end 
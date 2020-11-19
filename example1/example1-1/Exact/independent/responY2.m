function [output]=responY2(X,mu)
global cntY2
cntY2 = cntY2 + 1;
%multivariate Gaussian score function  
% shift transformation of X 
x1 = mu(1)+X(1);
x2 = mu(2)+X(2);
output = x1 + x2 - 6.45;
end 
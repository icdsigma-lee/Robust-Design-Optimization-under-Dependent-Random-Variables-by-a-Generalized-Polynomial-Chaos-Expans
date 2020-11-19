function [output]=responY1(X,mu)
%multivariate Gaussian score function  
% shift transformation of X 
global cntY1
cntY1 = cntY1 + 1;
x1 = mu(1)+X(1);
x2 = mu(2)+X(2);
output = (x1-4)^3 + (x1-3)^4 + (x2-5)^2+10;
end 
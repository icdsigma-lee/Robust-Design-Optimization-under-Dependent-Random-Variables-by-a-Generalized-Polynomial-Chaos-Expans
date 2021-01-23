%% ========================================================================
% Example 1.1: Function of constraint function and its grandient  
% Input: design variables (dv) Output: constraint value and its design sensitivities
% Analytic formulae for the first two moment and their sensitivities
% obtained by Mathematica 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [c, ceq, DC,  DCeq] = exactConfun(dv)
global cntCon 
cntCon = cntCon + 1;
sig = [0.4, 0.4];
cor = 0.4;
m1Y2 = -6.45+dv(1)+dv(2);
m2Y2 = 41.6025 + dv(1)^2 - ...
  12.9*dv(2) + dv(2)^2 + ...
  dv(1)*(-12.9 + 2.*dv(2)) + ...
  sig(1)^2 + 2.*cor*sig(1)*sig(2) + ...
  sig(2)^2;

meanY2 = m1Y2;
varY2 = m2Y2 - m1Y2^2;

sen1m = zeros(2,1);
sen1m(1) = 1;
sen1m(2) = 1;
sen2m = zeros(2,1);
sen2m(1) = -12.9 + 2*dv(1) + 2.*dv(2);
sen2m(2) = -12.9+2.*dv(1)+2*dv(2);

c = 3*sqrt(varY2)-meanY2;
ceq = [];
DC = ((3/2)*(1/sqrt(varY2))*(sen2m - 2*meanY2*sen1m)-sen1m);
DCeq = [];
disp(c)
end 

%% ========================================================================
% Example 1.2: Function of objective function and its grandient 
% Input: dv (design variables)
% The analytic formulae of the first two moments and their sensitivity
% obtained by Mathmatica 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [objValue, objGrad] = exactObjfun(dv)
global cntObj
cntObj = cntObj + 1;
sig = [dv(1)*0.15, dv(2)*0.15];
cor = 0;
m1Y1 = 52-11*dv(1)^3+dv(1)^4-10*dv(2)+dv(2)^2+42*sig(1)^2+3*sig(1)^4+6*dv(1)^2*(7+sig(1)^2)-3*dv(1)*(20+11*sig(1)^2)+sig(2)^2;
m2Y1 = 2704 - 22*dv(1)^7 + ...
dv(1)^8 - 20*dv(2)^3 + dv(2)^4 + ...
  7968*sig(1)^2 + 9564*sig(1)^4 + ...
  3075*sig(1)^6 + 105*sig(1)^8 + ...
  dv(1)^6*(205 + 28*sig(1)^2) - ...
  6*dv(1)^5*(174 + 77*sig(1)^2) + ... 
  1200*cor*sig(1)*sig(2) + ...
  660*cor*sig(1)^3*sig(2) + ...
  204*sig(2)^2 + 84*sig(1)^2* ...
   sig(2)^2 + 168*cor^2*sig(1)^2* ...
   sig(2)^2 + 6*sig(1)^4*sig(2)^2 + ...
  24*cor^2*sig(1)^4*sig(2)^2 + ...
  3*sig(2)^4 + 6*dv(2)^2*...
   (34 + 14*sig(1)^2 + sig(1)^4 + ...
    sig(2)^2) + dv(1)^4*...
   (3188 - 20*dv(2) + 2*dv(2)^2 + ...
    3075*sig(1)^2 + 210*sig(1)^4 + ...
    2*sig(2)^2) - 4*dv(2)*...
   (260 + 210*sig(1)^2 + ...
    15*sig(1)^4 + 60*cor*sig(1)*...
     sig(2) + 33*cor*sig(1)^3*sig(2) + ...
    15*sig(2)^2) - 2*dv(1)^3*...
   (3092 + 11*dv(2)^2 + ...
    5220*sig(1)^2 + 1155*sig(1)^4 + ...
    40*cor*sig(1)*sig(2) + ...
    11*sig(2)^2 - 2*dv(2)*...
     (55 + 4*cor*sig(1)*sig(2))) + ...
  3*dv(1)^2*(3075*sig(1)^4 + ...
    140*sig(1)^6 + 4*dv(2)^2*...
     (7 + sig(1)^2) + 220*cor*sig(1)*...
     sig(2) - 4*dv(2)*(70 + ...
      10*sig(1)^2 + 11*cor*sig(1)*...
       sig(2)) + ...
    4*(664 + 7*sig(2)^2) + ...
    4*sig(1)^2*(1594 + sig(2)^2 + ...
      2*cor^2*sig(2)^2)) - ...
  6*dv(1)*(2610*sig(1)^4 + ...
    385*sig(1)^6 + dv(2)^2*...
     (20 + 11*sig(1)^2) + ...
    280*cor*sig(1)*sig(2) + ...
    40*cor*sig(1)^3*sig(2) - ...
    2*dv(2)*(100 + 55*sig(1)^2 + ...
      28*cor*sig(1)*sig(2) + ...
      4*cor*sig(1)^3*sig(2)) + ...
    20*(52 + sig(2)^2) + ...
    sig(1)^2*(3092 + ...
      11*(1 + 2*cor^2)*sig(2)^2));
meanY1 = m1Y1;
varY1 = m2Y1 - m1Y1^2;

sen1m = zeros(2,1);
sen1m(1) = -60 - ...
  33*dv(1)^2 + 4*dv(1)^3 - ...
  33*sig(1)^2 + 12*dv(1)*...
   (7 + sig(1)^2);
sen1m(2) = 2*(-5 + dv(2));
sen2m = zeros(2,1);
sen2m(1) = -154*dv(1)^6 + 8*dv(1)^7 + ...
  6*dv(1)^5*(205 + 28*sig(1)^2) - ...
  30*dv(1)^4*(174 + 77*sig(1)^2) + ...
  4*dv(1)^3*(3188 - 20*dv(2) + ...
    2*dv(2)^2 + 3075*sig(1)^2 + ...
    210*sig(1)^4 + 2*sig(2)^2) - ...
  6*dv(1)^2*(3092 + 11*dv(2)^2 +... 
    5220*sig(1)^2 + 1155*sig(1)^4 +... 
    40*cor*sig(1)*sig(2) + ...
    11*sig(2)^2 - 2*dv(2)*...
     (55 + 4*cor*sig(1)*sig(2))) +... 
  6*dv(1)*(3075*sig(1)^4 + ...
    140*sig(1)^6 + 4*dv(2)^2*...
     (7 + sig(1)^2) + 220*cor*sig(1)*...
     sig(2) - 4*dv(2)*(70 + ...
      10*sig(1)^2 + 11*cor*sig(1)*...
       sig(2)) + ...
    4*(664 + 7*sig(2)^2) + ...
    4*sig(1)^2*(1594 + sig(2)^2 + ...
      2*cor^2*sig(2)^2)) - ...
  6*(2610*sig(1)^4 + 385*sig(1)^6 + ...
    dv(2)^2*(20 + 11*sig(1)^2) + ...
    280*cor*sig(1)*sig(2) + ...
    40*cor*sig(1)^3*sig(2) - ...
    2*dv(2)*(100 + 55*sig(1)^2 + ...
      28*cor*sig(1)*sig(2) + ...
      4*cor*sig(1)^3*sig(2)) + ...
    20*(52 + sig(2)^2) + ...
    sig(1)^2*(3092 + ...
      11*(1 + 2*cor^2)*sig(2)^2));
sen2m(2) = 4* dv(1)^4*(-5 + dv(2)) - ...
  60*dv(2)^2 + 4*dv(2)^3 - ...
  4*dv(1)^3*(-55 + 11*dv(2) - ...
    4*cor*sig(1)*sig(2)) + ...
  12*dv(1)*(100 + 55*sig(1)^2 - ...
    dv(2)*(20 + 11*sig(1)^2) + ...
    28*cor*sig(1)*sig(2) + ...
    4*cor*sig(1)^3*sig(2)) + ...
  12*dv(2)*(34 + 14*sig(1)^2 + ...
    sig(1)^4 + sig(2)^2) - ...
  4*(260 + 210*sig(1)^2 + ...
    15*sig(1)^4 + 60*cor*sig(1)*...
     sig(2) + 33*cor*sig(1)^3*sig(2) + ...
    15*sig(2)^2) + 3*dv(1)^2*...
   (8*dv(2)*(7 + sig(1)^2) - ...
    4*(70 + 10*sig(1)^2 + ...
      11*cor*sig(1)*sig(2)));
  objValue = sqrt(varY1)/45;
  objGrad = ((1/90)*(1/sqrt(varY1))*(sen2m-2*meanY1*sen1m)); 
  disp(varY1)
end 

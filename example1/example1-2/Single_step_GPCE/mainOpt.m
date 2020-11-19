%% ========================================================================
%  main program (example1-2)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% RDO for mathematical functions
% save history for every iterations 
% data name 
% - Exact : resultE
% - Direct approach : resultD
% - Singular Step GPCE : resultS
% - Multipoint approximation : resultM
%% ========================================================================
clear all
clc 
global cntObj cntCon cntY1 cntY2
cntObj=0; cntCon=0; cntY1=0; cntY2=0;
x0 = [5,5];
estS = cell(5,1);

[stat0] = estObjfun(x0);
[c0,~,~,~] =  confun(x0);
FilNam = sprintf('resultS41.mat');
[historyS, searchdirS] = runfmincon(x0);
estS{5,1} = [cntY1, cntY2]; % function call # for Y1 and Y2  

xf = historyS.x(end,:);
cntObj=0;
cntCon=0;
[statf] = estObjfun(xf);
[cf,~,~,~] =  confun(xf);

estS{1,1}  = stat0; % mean and variance at initial design 
estS{2,1} = statf; % mean and variance at optimum design 
estS{3,1} = c0; % constraint value at initial design
estS{4,1}  = cf; % constraint value at optimum design 

save(FilNam, 'historyS', 'estS');
%iter=1:length(history.fval);
%plot(iter,history.fval,);
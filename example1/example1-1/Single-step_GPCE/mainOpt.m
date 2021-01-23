%% ========================================================================
%  main program (example1-1)
% RDO for mathematical functions
% save history for every iterations 
% output data title 
% - Exact : resultE
% - Direct approach : resultD
% - Singular Step GPCE : resultS
% - Multipoint approximation : resultM
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
clear all
clc 
global cntObj cntCon 
global cntY1 cntY2
global stat0 statf sopt

estS = cell(5,1);
x0 = [5,5]; % initial design 

%% monomial moment matrix (Gm)
genGramMatrix % Save ID (multi-index j), QQ (whitenning tran. matrix)

cntObj=0; cntCon=0; cntY1 = 0; cntY2 = 0;

%% options: 'pre', 'run', 'post' 
% 'pre' : save the state of initial design 
% 'run' : run optimization  
% 'post' : save the state of optimal design 

% estimation of the state of objective/constraint functions at initial design
sopt = 'pre';
[~] = objfun(x0);
[c0,~,~,~] =  confun(x0);

FilNam = sprintf('resultS.mat');

% run optimization
sopt = 'run';
[historyS, searchdirS] = runfmincon(x0);
estS{5,1} = [cntY1, cntY2]; % # of function evaluations for Y1 and Y2  
xf = historyS.x(end,:);

% estimation of the state of objective/constraint functions at optimal design
sopt = 'post';  

[~,~] = objfun(xf); 
[cf,~,~,~] =  confun(xf);

estS{1,1}  = stat0; % mean and variance at initial design 
estS{2,1} = statf; % mean and variance at optimal design 
estS{3,1} = c0; % constraint functions at initial design
estS{4,1}  = cf; % constraint functions at optimum design 
save(FilNam, 'historyS','estS');
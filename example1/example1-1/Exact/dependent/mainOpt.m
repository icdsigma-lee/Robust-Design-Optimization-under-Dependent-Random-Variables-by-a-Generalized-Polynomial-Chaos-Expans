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

estExact = cell(5,1);
x0 = [5,5];

cntObj=0; cntCon=0; 

%% options: 'pre', 'run', 'post' 
% 'pre' : save the state of initial design 
% 'run' : run optimization  
% 'post' : save the state of optimal design 

% estimation of the state of objective/constraint functions at initial design 
sopt = 'pre';
[stat0] = estExactobjfun(x0);
[c0,~,~,~] =  exactConfun(x0);

FilNam = sprintf('resultExact.mat');
cntY1 = 0; cntY2 = 0;

% run optimization
sopt = 'run';
[historyExact, searchdirExact] = runfmincon(x0);
estExact{5,1} = [cntY1, cntY2]; % # of function evaluations for Y1 and Y2  
xf = historyExact.x(end,:);

% estimation of the state of objective/constraint functions at optimal design  
sopt = 'post';  

[statf] = estExactobjfun(xf); 
[cf,~,~,~] =  exactConfun(xf);

estExact{1,1}  = stat0; % mean and variance at initial design 
estExact{2,1} = statf; % mean and variance at optimal design 
estExact{3,1} = c0; % constraint functions at initial design
estExact{4,1}  = cf; % constraint functions at optimal design 
save(FilNam, 'historyExact','estExact');

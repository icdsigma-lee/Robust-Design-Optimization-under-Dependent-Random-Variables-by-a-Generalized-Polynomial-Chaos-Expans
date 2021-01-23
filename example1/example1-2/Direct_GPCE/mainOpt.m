%% ========================================================================
%  main program (example1-2)
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

cntObj=0; cntCon=0; cntY1=0; cntY2=0;
estD = cell(5,1);
x0 = [5,5]; % initial design 

%% options: 'pre', 'run', 'post' 
% 'pre' : save the state of initial design 
% 'run' : run optimization  
% 'post' : save the state of optimal design 

% estimation of the state of objective/constraint functions at initial design
sopt = 'pre';
[stat0] = estObjfun(x0);
[c0,~,~,~] =  confun(x0);

FilNam = sprintf('resultD.mat');
cntY1 = 0; cntY2 = 0;

% run optimization
sopt = 'run';
[historyD, searchdirD] = runfmincon(x0);
estD{5,1} = [cntY1, cntY2]; % # of function evaluations for Y1 and Y2  
xf = historyD.x(end,:);

% estimation of the state of objective/constraint functions at optimal design 
sopt = 'post';  
[statf] = estObjfun(xf);
[cf,~,~,~] =  confun(xf);

estD{1,1}  = stat0; % mean and variance at initial design 
estD{2,1} = statf; % mean and variance at optimal design 
estD{3,1} = c0; % constraint value at initial design
estD{4,1}  = cf; % constraint value at optimum design 
save(FilNam, 'historyD','estD');

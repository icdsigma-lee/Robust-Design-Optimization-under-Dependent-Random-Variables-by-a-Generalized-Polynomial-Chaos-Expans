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
global cntObj cntCon 
global cntY1 cntY2
global stat0 statf sopt

cntObj=0; cntCon=0; cntY1=0; cntY2=0;
estD = cell(5,1);

x0 = [5,5]; % initial design 

%% objective has options as 'pre', 'run', 'post' 
% 'pre' : save the stat# at the initial design 
% 'run' : run optimization  
% 'post' : save the stat# at the optimum design

% pre-estimation of constraint value at init design
sopt = 'pre';
[stat0] = estObjfun(x0);
[c0,~,~,~] =  confun(x0);
FilNam = sprintf('resultD.mat');
cntY1 = 0; cntY2 = 0;
% run optimization
sopt = 'run';
[historyD, searchdirD] = runfmincon(x0);
estD{5,1} = [cntY1, cntY2]; % function call # for Y1 and Y2  
xf = historyD.x(end,:);
% post-estimation of statf and constraint value for optimum  
sopt = 'post';  

[statf] = estObjfun(xf);
[cf,~,~,~] =  confun(xf);

estD{1,1}  = stat0; % mean and variance at initial design 
estD{2,1} = statf; % mean and variance at optimum design 
estD{3,1} = c0; % constraint value at initial design
estD{4,1}  = cf; % constraint value at optimum design 
save(FilNam, 'historyD','estD');

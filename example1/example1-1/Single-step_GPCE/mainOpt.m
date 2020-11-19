%% ========================================================================
%  main program (example1-1)
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

estS = cell(5,1);
x0 = [5,5];

%% create an monomial moment matrix for initial design
genGramMatrix % Save ID, QQ (whitenning tran. matrix)

cntObj=0; cntCon=0; cntY1 = 0; cntY2 = 0;

%% objective has options as 'pre', 'run', 'post' 
% 'pre' : save the stat# at the initial design 
% 'run' : run optimization  
% 'post' : save the stat# at the optimum design

% pre-estimation of constraint value for init. 
sopt = 'pre';
[~] = objfun(x0);
[c0,~,~,~] =  confun(x0);

FilNam = sprintf('resultS.mat');
% run optimization
sopt = 'run';
[historyS, searchdirS] = runfmincon(x0);
estS{5,1} = [cntY1, cntY2]; % function call # for Y1 and Y2  
xf = historyS.x(end,:);
% post-estimation of stat and constraint value for optimum  
sopt = 'post';  

[~,~] = objfun(xf); 
[cf,~,~,~] =  confun(xf);

estS{1,1}  = stat0; % mean and variance at initial design 
estS{2,1} = statf; % mean and variance at optimum design 
estS{3,1} = c0; % constraint value at initial design
estS{4,1}  = cf; % constraint value at optimum design 
save(FilNam, 'historyS','estS');
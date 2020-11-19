%% ========================================================================
%  main program (example3)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% RDO for mathematical functions
% save history for every iterations 
% data name 
% - Exact : resultE (Finite difference method)
% - Direct approach : resultD# (#th-order GPCE)
% - Singular Step GPCE : resultS# (#th-order GPCE)
% - Multipoint approximation : resultM# (#th-order GPCE)
%% ========================================================================
clear all
clc 

global cntObj cntCon  
global cntY0 cntY1 cntY2
global stat0 statf
global sopt 

estE = cell(6,1);
x0 = [20,20,1,1]; % initial design 

%% estimate statistics and do feasible estimation  
% Obj. and const. func. saves the relevant info. at the first iter.  
cntObj=0; cntCon=0; cntY0=0; cntY1=0; cntY2=0;

%% objective has options as 'pre', 'run', 'post' 
% 'pre' : save the stat. for the initial design 
% 'run' : do not save any stat. info. 
% 'post' : save the stat. for the optimum 

% pre-estimation of stat and constraint value for init. 
sopt = 'pre';
[~,~] = objfun(x0); 
[c0,~,~,~] =  confun(x0);

% init. count # 
cntObj=1; cntCon=1; cntY0=0; cntY1=0; cntY2=0;
FilNam = sprintf('resultE.mat');
% run optimization
sopt = 'run';
tic;

[historyE, searchdir] = runfmincon(x0);
tE = toc;
xf = historyE.x(end,:);
estE{5,1} = [cntY0 cntY1, cntY2]; % function call # for Y1 and Y2 
% post-estimation of stat and constraint value for optimum  
sopt = 'post';
[~,~] = objfun(xf); 
[cf,~,~,~] =  confun(xf);
estE{1,1}  = stat0; % mean and variance at initial design 
estE{2,1} = statf; % mean and variance at optimum design 
estE{3,1} = c0; % constraint value at initial design
estE{4,1}  = cf; % constraint value at optimum design 
estE{6,1}  = tE;% CPU time 
save(FilNam, 'historyE','estE');


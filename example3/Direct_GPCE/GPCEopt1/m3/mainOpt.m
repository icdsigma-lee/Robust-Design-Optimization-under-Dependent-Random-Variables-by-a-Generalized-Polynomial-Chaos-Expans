%% ========================================================================
%  main program (example2)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% RDO for mathematical functions
% save history for every iterations 
% data name 
% - Exact : resultE
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

estD3 = cell(6,1);
x0 = [20,20,1,1]; % initial design 

%% create an monomial moment matrix for initial design
%genGramMatrix; % Save ID, G1, QQ (whitenning tran. matrix)
%modGramMatrix; % Save ID, INDEX#, QQ# (sub whitenning tran. matrix)  

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
FilNam = sprintf('resultD3.mat');
% run optimization
sopt = 'run';
tic;
[historyD3, searchdirD3] = runfmincon(x0);
tD3 = toc;
xf = historyD3.x(end,:);
estD3{5,1} = [cntY0 cntY1, cntY2]; % function call # for Y1 and Y2  
% post-estimation of stat and constraint value for optimum  
sopt = 'post';
[~,~] = objfun(xf); 
[cf,~,~,~] =  confun(xf);
estD3{1,1}  = stat0; % mean and variance at initial design 
estD3{2,1} = statf; % mean and variance at optimum design 
estD3{3,1} = c0; % constraint value at initial design
estD3{4,1}  = cf; % constraint value at optimum design 
estD3{6,1}  = tD3;% CPU time 
save(FilNam, 'historyD3','estD3');
%% ========================================================================
%  main program (example4) 
% RDO problem
% save history for every iterations 
% data name 
% - Exact : resultE
% - Direct approach : resultD# (#th-order GPCE)
% - Singular Step GPCE : resultS# (#th-order GPCE)
% - Multipoint approximation : resultM# (#th-order GPCE)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu)
%% ========================================================================
clear all
clc 

global cntObj cntCon  
global cntRspObj cntRspCon
global stat0 statf
global sopt 

estD2 = cell(5,1);
x0 = ones(1,10)*30; %initial design 

%% Monomial moment matrix 
genGramMatrix % Save ID, QQ (whitenning tran. matrix)

cntObj=0; cntCon=0; cntRspObj=0; cntRspCon=0;

%% options: 'pre', 'run', 'post' 
% 'pre' : save the state of initial design 
% 'run' : run optimization  
% 'post' : save the state of optimal design 

% estimation of the state of objective/constraint functions at initial design 
sopt = 'pre';
[~,~] = objfun(x0); 
[c0,~,~,~] =  confun(x0);

% init. count # 
FilNam = sprintf('resultD2.mat');
% run optimization
sopt = 'run';
[historyD2, searchdirD2] = runfmincon(x0);
xf = historyD2.x(end,:);
estD2{5,1} = [ cntRspObj, cntRspCon]; % function call # for Y1 and Y2  
% estimation of the state of objective/constraint functions at optimal design  
sopt = 'post';
[~,~] = objfun(xf); 
[cf,~,~,~] =  confun(xf);
estD2{1,1}  = stat0; % mean and variance at initial design 
estD2{2,1} = statf; % mean and variance at optimal design
estD2{3,1} = c0; % constraint value at initial design
estD2{4,1}  = cf; % constraint value at optimal design 
save(FilNam, 'historyD2','estD2');
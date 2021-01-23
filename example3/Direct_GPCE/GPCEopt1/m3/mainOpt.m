%% ========================================================================
%  main program (example3)
% RDO for mathematical functions
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
global cntY0 cntY1 cntY2
global stat0 statf
global sopt 

estD3 = cell(6,1);
x0 = [20,20,1,1]; % initial design 

%% Monomial moment matrix 
genGramMatrix; % Save ID, G1, QQ (whitening tran. matrix)
modGramMatrix; % Save ID, INDEX#, QQ# (sub whitening tran. matrix)  

   
cntObj=0; cntCon=0; cntY0=0; cntY1=0; cntY2=0;

%% options: 'pre', 'run', 'post' 
% 'pre' : save the state of initial design 
% 'run' : run optimization  
% 'post' : save the state of optimal design 

% estimation of the state of objective/constraint functions at initial design 
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
% estimation of the state of objective/constraint functions at optimal design  
sopt = 'post';
[~,~] = objfun(xf); 
[cf,~,~,~] =  confun(xf);
estD3{1,1}  = stat0; % mean and variance at initial design 
estD3{2,1} = statf; % mean and variance at optimal design   
estD3{3,1} = c0; % constraint value at initial design
estD3{4,1}  = cf; % constraint value at optimal design  
estD3{6,1}  = tD3;% CPU time 
save(FilNam, 'historyD3','estD3');
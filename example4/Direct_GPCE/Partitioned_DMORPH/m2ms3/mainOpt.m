%% ========================================================================
%  main program (example4)
% RDO problem
% save history for every iterations 
% data name 
% - Exact : resultE
% - Direct approach w/SLS : resultD# (#th-order GPCE)
% - Direct approach w/partitioned DMORPH : resultDDM# (#th-order GPCE)
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

est = cell(5,1);
x0 = ones(1,10)*30; % initial design 

%% Monomial moment matrix 
genGramMatrix % Save ID, QQ (whitening tran. matrix)

 
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
FilNam = sprintf('resultDDM2.mat');
% run optimization
sopt = 'run';
[history, searchdir] = runfmincon(x0);
xf = history.x(end,:);
cntRspObj1 = cntRspObj;
cntRspCon1 = cntRspCon;
est{5,1} = [ cntRspObj1, cntRspCon1]'; % function call # for Y1 and Y2  
% estimation of the state of objective/constraint functions at optimal design  
sopt = 'post';
[~,~] = objfun(xf); 
[cf,~,~,~] =  confun(xf);
est{1,1}  = stat0'; % mean and variance at initial design 
est{2,1} = statf'; % mean and variance at optimal design  
est{3,1} = c0'; % constraint value at initial design
est{4,1}  = cf'; % constraint value at optimal design  
save(FilNam, 'history','est');
%% ========================================================================
%  main program (example3)
% RDO for mathematical functions
% save history for every iterations 
% data name 
% - Exact : resultE
% - Direct approach (opt1) : resultD# (#th-order GPCE)
% - Direct approach (opt2) : resultD#opt2 (#th-order GPCE)
% - Singular Step GPCE : resultS# (#th-order GPCE)
% - Multipoint approximation : resultM# (#th-order GPCE)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
clear all
clc 
global cntObj cntCon  
global cntY0 cntY1 cntY2
global stat0 statf
global sopt ii jj initobj initcon

est = cell(6,1);
x0 = [20,20,1,1]; %initial design 

%% Monomial moment matrix 
genGramMatrix; % Save ID, QQ (whitenning tran. matrix)
modGramMatrix;
  
cntObj=0; cntCon=0; cntY0=0; cntY1=0; cntY2=0;
ii = 0; jj = 0;
initobj = 0;
initcon = 0;

%% options: 'pre', 'run', 'post' 
% 'pre' : save the state of initial design 
% 'run' : run optimization  
% 'post' : save the state of optimal design 

% estimation of the state of objective/constraint functions at initial design
sopt = 'pre';
[~,~] = sobjfun(x0); 
[c0,~,~,~] =  sconfun(x0);
cntY0=0; cntY1=0; cntY2=0;
% init. count # 
FilNam = sprintf('resultM1.mat');
% run optimization
sopt = 'run';
tic;
[history, searchdir] = runfmincon(x0);
tM1 = toc; 
xf = history.x(end,:)';
est{5,1} = [cntY0, cntY1, cntY2]; % function call # for Y1 and Y2  
% estimation of the state of objective/constraint functions at optimal design  
sopt = 'post';
[~,~] = sobjfun(xf); 
[cf,~,~,~] =  sconfun(xf);
est{1,1}  = stat0'; % mean and variance at initial design 
est{2,1} = statf'; % mean and variance at optimal design 
est{3,1} = c0'; % constraint value at initial design
est{4,1}  = cf'; % constraint value at optimal design 
est{6,1}  = tM1;
save(FilNam, 'history','est');
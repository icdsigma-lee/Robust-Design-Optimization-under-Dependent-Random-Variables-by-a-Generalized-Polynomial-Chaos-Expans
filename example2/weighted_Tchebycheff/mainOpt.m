%% ========================================================================
%  main program (example2)
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
global stat0 statf sopt w1 w2

estS = cell(5,1);
x0 = [5,5];

%% monomial moment matrix (Gm) 
genGramMatrix  % Save ID (multi-index j), QQ (whitenning tran. matrix)

cntObj=0; cntCon=0; cntY1 = 0; cntY2 = 0;
vw1 = [0.1:0.1:0.9];
vw2 = [0.9:-0.1:0.1];
lenVw = length(vw1);
for i=1:lenVw 
%% options: 'pre', 'run', 'post' 
% 'pre' : save the state of initial design 
% 'run' : run optimization  
% 'post' : save the state of optimal design 
w1 = vw1(i);
w2 = vw2(i);
cntY1 = 0; cntY2 = 0;
% estimation of the state of objective/constraint functions at initial design
sopt = 'pre';
[~] = objfun(x0);
[c0,~,~,~] =  confun(x0);

FilNam = sprintf('resultS%f%f.mat',w1,w2);
% run optimization
sopt = 'run';
[historyS, searchdirS] = runfmincon(x0);
estS{5,1} = [cntY1, cntY2]; % # of function evaluations for Y1 and Y2  
xf = historyS.x(end,:);

% estimation of the state of objective/constraint functions at optimal design   
sopt = 'post';  
[~,~] = objfun(xf); 
[cf,~,~,~] =  confun(xf);

estS{1,1}  = stat0; % mean and variance at initial design 
estS{2,1} = statf; % mean and variance at optimal design 
estS{3,1} = c0; % constraint functions at initial design
estS{4,1}  = cf; % constraint functions at optimal design 
save(FilNam, 'historyS','estS');
end 
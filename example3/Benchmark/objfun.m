%% ========================================================================
% Example 2: Function of objective function and its grandient (centrial finite difference method)   
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Input required: design variables (dv) 
%% ========================================================================
function [objValue, objGrad] = objfun(dv)
global cntObj stat0 statf sopt

double precision;
%% Initialization
% number of variables
cntObj = cntObj + 1;
FilNam = sprintf('data.mat');
% Define parameters of distributiosn
N = 7;
nd = 4;
mu5 = 10000; %kgm^3
mu6 = 2050000000; %Pa
mu7 = 200000; %N
% normalized mean vector for N variables 
mu = ones(1,N);
% transform to original mean vector 
muTr = [dv(1), dv(2), dv(3), dv(4), mu5, mu6, mu7];
diff = 1E-6;
% standard deviatiaon from
sig = [0.02, 0.02, 0.02, 0.02, 0.3, 25/105, 0.25];

% Correlation coefficient matrix
cor = zeros(N,N);
for i=1:N
    for j=i:N
        if ((i==1) && (j==2))
            cor(i,j) = 0.4;
        end 
        if ((i==3) && (j==4))
            cor(i,j) = -0.4;
        end 
        if (i==j)
            cor(i,j) = 1;
        end 
    end 
end 
for i=1:N
    for j=i:N
        cor(j,i) = cor(i,j);
    end 
end 
nSample0 = 5000000;
nSample = 500000;
if (cntObj == 1)
%% Sampling

% Sample generation for X1~X7
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z = qrand(q,nSample0);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);

% transformation for x1~x2
cov12 = zeros(2,2);
for i = 1:2
    for j = i:2
        cov12(i,j) = cor(i,j)*sig(i)*sig(j);
    end 
end 
for i = 1:2
    for j=i:2
        cov12(j,i) = cov12(i,j);
    end 
end 
T12 = chol(cov12,'lower');
x(:,1:2) = (T12*x(:,1:2)')';
for i=1:2
    x(:,i) = x(:,i) + mu(i);
end 
% transformation for x3~x4
cov34 = zeros(2,2); %match index by + 2 
for i = 1:2
    for j = i:2
        cov34(i,j) = cor(i+2,j+2)*sig(i+2)*sig(j+2);
    end 
end 
for i = 1:2
    for j=i:2
        cov34(j,i) = cov34(i,j);
    end 
end 
T34 = chol(cov34,'lower');
x(:,3:4) = (T34*x(:,3:4)')';
for i=1:2
    x(:,i+2) = x(:,i+2) + mu(i+2);
end 
cov = zeros(nd,nd);
cov(1:2,1:2) = cov12;
cov(3:4,3:4) = cov34;
%% transformation for x5 (weibull distribution) 
k = 3.71377203781000;
lamda = 1.10786387179285;
x(:,5) = (lamda)*(-log(1-normcdf(x(:,5)))).^(1/k);
%% transformation for x6 (gumbel distribution) 
beta6 = 0.185642095531828;
mu6 = 0.89284447439388222364137325872783;
x(:,6) = mu6 + beta6*(-log(-log(normcdf(x(:,6)))));
%% transformation for x7 (gumbel distribution)
beta7 = 0.194924200308419;
mu7 = 0.88748669811357634764020434862424;
x(:,7) = mu7 + beta7*(-log(-log(normcdf(x(:,7)))));

tmpY = zeros(nSample,1);
% generate data 
x(nSample+1:end,:) = [];
	for L = 1:nSample 
        tmpY(L,1) = responY0(x(L,:),muTr);
	end 
%% Sensitivity analysis
% first moment
sen1m = zeros(1,nd);
% second moment 
sen2m = zeros(1,nd);
tmpYL = zeros(nSample,1);
tmpYU = zeros(nSample,1);
for i = 1:nd
 	% reset muUpdateL and muUpdateU
  		muUpdateL = muTr;
  		muUpdateU = muTr;
    if (i < 3)
    	muUpdateL(i) = dv(i)-diff; 
    	muUpdateU(i) = dv(i)+diff;
    else 
    	muUpdateL(i) = dv(i) - diff;
    	muUpdateU(i) = dv(i) + diff;
    end
    
    for L=1:nSample 
    	tmpYL(L,1) = responY0(x(L,:), muUpdateL);
    	tmpYU(L,1) = responY0(x(L,:), muUpdateU);   
    end     
    sen1m(1,i) = (mean(tmpYU) - mean(tmpYL))/(2*diff);
    sen2m(1,i) = (mean(tmpYU.^2) - mean(tmpYL.^2))/(2*diff);
end 

varY0 = var(tmpY);
meanY0 = mean(tmpY);

w1 = meanY0;
w2 = sqrt(varY0);

objValue = 0.5*meanY0/w1 + 0.5*sqrt(varY0)/w2;
objGrad = ((0.5/w1)*sen1m + (0.5/w2)*(1/sqrt(varY0))*(sen2m - 2*meanY0*sen1m))'; 
disp(dv) 
disp(varY0)
disp(meanY0)
save(FilNam, 'x','w1','w2');



else % cntObj ~= 1
load(FilNam);
tmpY = zeros(nSample,1);
for L = 1:nSample 
        tmpY(L,1) = responY0(x(L,:),muTr);
end 
%% Sensitivity analysis
% first moment
sen1m = zeros(1,nd);
% second moment 
sen2m = zeros(1,nd);
tmpYL = zeros(nSample,1);
tmpYU = zeros(nSample,1);
for i = 1:nd
 	% reset muUpdate 
  		muUpdateL = muTr;
  		muUpdateU = muTr;
    if (i < 3)
    	muUpdateL(i) = dv(i)-diff; 
    	muUpdateU(i) = dv(i)+diff;
    else 
    	muUpdateL(i) = dv(i) - diff;
    	muUpdateU(i) = dv(i) + diff;
    end
    
    for L=1:nSample 
    	tmpYL(L,1) = responY0(x(L,:), muUpdateL);
    	tmpYU(L,1) = responY0(x(L,:), muUpdateU);
    end 
        
    sen1m(1,i) = (mean(tmpYU) - mean(tmpYL))/(2*diff);
    sen2m(1,i) = (mean(tmpYU.^2) - mean(tmpYL.^2))/(2*diff);
end 

varY0 = var(tmpY);
meanY0 = mean(tmpY);

objValue = 0.5*meanY0/w1 + 0.5*sqrt(varY0)/w2;
objGrad = ((0.5/w1)*sen1m + (0.5/w2)*(1/sqrt(varY0))*(sen2m - 2*meanY0*sen1m))'; 
disp(dv) 
disp(varY0)
disp(meanY0)

end % end of function 

% record stat. info. at the initial and final design. 
switch sopt
    case 'pre'
        stat0 = [meanY0, sqrt(varY0)];
    case 'post'
        statf = [meanY0, sqrt(varY0)];
    otherwise 
end 




            

















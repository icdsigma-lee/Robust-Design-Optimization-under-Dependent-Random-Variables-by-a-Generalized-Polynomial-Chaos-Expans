%% ========================================================================
% Example 3: Function of inequality constraint function and its grandient (central finite difference method)   
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Input required: design variables (dv) 
%% ========================================================================
function[c, ceq, DC,  DCeq] = confun(dv)

global cntCon con
double precision;
%% Initialization
% number of variables
cntCon = cntCon + 1;
FilNam = sprintf('data.mat');
% Define parameters of distributiosn 
N = 7;
nd  = 4;
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
nSample = 500000;
load(FilNam); 
tmpY1 = zeros(nSample,1);
tmpY2 = zeros(nSample,1);

for L=1:nSample
    	tmpY1(L,1) = responY1(x(L,:),muTr);
        tmpY2(L,1) = responY2(x(L,:),muTr);
end 
 
%% Sensitivity analysis 
% first moment
sen1m1 = zeros(1,nd);
sen1m2 = zeros(1,nd);
sen2m1 = zeros(1,nd);
sen2m2 = zeros(1,nd);
tmpY1L = zeros(nSample,1);
tmpY1U = zeros(nSample,1);
tmpY2L = zeros(nSample,1);
tmpY2U = zeros(nSample,1);

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
    	tmpY1L(L,1) = responY1(x(L,:), muUpdateL);
    	tmpY1U(L,1) = responY1(x(L,:), muUpdateU);
    	tmpY2L(L,1) = responY2(x(L,:), muUpdateL);
    	tmpY2U(L,1) = responY2(x(L,:), muUpdateU);
    end 
        
    sen1m1(1,i) = (mean(tmpY1U) - mean(tmpY1L))/(2*diff);
    sen1m2(1,i) = (mean(tmpY2U) - mean(tmpY2L))/(2*diff);
    sen2m1(1,i) = (mean(tmpY1U.^2) - mean(tmpY1L.^2))/(2*diff);
    sen2m2(1,i) = (mean(tmpY2U.^2) - mean(tmpY2L.^2))/(2*diff);
end 

varY1 = var(tmpY1);
meanY1 = mean(tmpY1);
varY2 = var(tmpY2);
meanY2 = mean(tmpY2);
c(1) = meanY1 + 3*sqrt(varY1);
c(2) = meanY2 + 3*sqrt(varY2);
ceq = [];
DC1 = sen1m1' + (3/2)*(1/sqrt(varY1))*(sen2m1-2*meanY1*sen1m1)';
DC2 = sen1m2' + (3/2)*(1/sqrt(varY2))*(sen2m2-2*meanY2*sen1m2)';
DC = [DC1 , DC2];
DCeq = [];
disp('inequality:')
disp(c)
con = c; % save c to global con 
% End of function 







            

















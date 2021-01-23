%% ========================================================================
% Example 1.1: Function of inequality constraint function and its grandient  
% Input: design variables (dv) 
% Output: constraint value and its design sensitivities
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [c, ceq, DC,  DCeq] = confun(dv)
global cntCon
double precision;
%% Initialization
cntCon = cntCon + 1; % count function call #
N = 2; % # of random variables 
m = 1;  % order of GPCE for generic function 
ms = 1; % order of GPCE for score function
nd = 2; % # of design variables 
nA = nchoosek(N+m, m); % # of GPCE expansion coefficients for generic function
nAs = nchoosek(N+ms, ms);% # of GPCE expansion coefficients for score function
grA = [nA, nAs];
A = max(grA);
% files for data at 1st iteration 
FilNam = sprintf('data.mat');
FilNam1 = sprintf('data1.mat');
% means vector (mu) 
mu1 = 1; 
mu2 = 1;
mu = [mu1, mu2];
% standard deviation (sig) 
sig1 = 0.15;
sig2 = 0.15;
sig = [sig1, sig2];
% correlation matrix  
rho12 = -0.5;
cor = zeros(N,N);
% covariance matrix (cov)
cov = zeros(N,N);
for i=1:N
    for j=1:N
        if (i==j)
            cov(i,i) = sig(i)^2;
        else 
            cov(i,j) = rho12*sig(i)*sig(j);
        end 
            cor(i,j) = cov(i,j)/(sig(i)*sig(j));
    end 
end 
% original covariance matrix (covT) 
covT = zeros(N,N);
for i=1:N
    for j=1:N
        if (i==j)
            covT(i,i) = (dv(i)*sig(i))^2;
        else 
            covT(i,j) = rho12*dv(i)*sig(i)*dv(j)*sig(j);
        end 
    end 
end 

nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 1000000;
grSample = [nSampley, nSamples, nSampleo];
nSample = max(grSample); 
% first iteration 
if (cntCon == 1)
load(FilNam);   
% output data 
for L = 1:nSample 
    if (L < nSampley +1) 
        tmpY(L,1) = responY2(x(L,:),dv);
    end 
end 

infoMy = infoM([1:nSampley],[1:nA]);
infoMs = infoM([1:nSamples],[1:nAs]);
cs = zeros(nAs,nd);
cy = (infoMy'*infoMy)\(infoMy'*tmpY);
for i=1:nd 
    cs(:,i) = (infoMs'*infoMs)\(infoMs'*tmpS(:,i));
end 
cs(1,:) = 0;
%% Sensitivity analysis 
% first moment
 jacob = zeros(nd,nd);
for i = 1:nd
    jacob(i,i) = 1/dv(i);
end 

sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs(j,i);
    end
end 
sen1m = sen1m*jacob;
% second moment 
sen2m = zeros(1,nd);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1))
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);     end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);     end 
                end 
            else 
                % E[ON X ON X ON]          
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*triON(i,j,k);
            end 
        end 
    end 
 end  
sen2m = sen2m*jacob;
% variance 
varY2 = sum(cy(2:end).^2);
% mean
meanY2 = cy(1);
% constraint value 
c = 3*sqrt(varY2)-meanY2;
ceq = [];
% design sensitivities
DC = ((3/2)*(1/sqrt(varY2))*(sen2m - 2*meanY2*sen1m)-sen1m)';
DCeq = [];
cy1 = cy;
save(FilNam1, 'cy1','cs');
% the rest of the first iteration
% single-step GPCE 
else 
load(FilNam);
load(FilNam1);
zTran = zeros(nd, nd);
for i = 1:nd
    zTran(i,i) = dv(i)/dv0(i);
end 
x(nSample+1:end,:) = [];
x = x*zTran;
% monomial basis (M) 
M = zeros(nSample,A); 
for i=1:A
    chkID = ID(i,:);
    nZeroID = find(chkID~=0); 
    nZero = length(nZeroID);
    if (nZero == 0)
        M(:,i) = 1;
    end 
    if (nZero == 1)
        M(:,i) = x(:,nZeroID).^chkID(nZeroID);
    end 
    if (nZero == 2)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        M(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2));
    end 
end  

% new output data using GPCE expansion coefficients at old design 
for L = 1:nSample 
    if (L < nSampley +1) 
        tmpY(L,1) = cy1'*QQ(1:nA,1:nA)*M(L,:)';
    end 
end 

infoMy = infoM([1:nSampley],[1:nA]);
infoMs = infoM([1:nSamples],[1:nAs]);

cy = (infoMy'*infoMy)\(infoMy'*tmpY);

%% Sensitivity analysis 
% first moment
 jacob = zeros(nd,nd);
for i = 1:nd
    jacob(i,i) = 1/dv(i);
end 

sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs(j,i);
    end
end 
sen1m = sen1m*jacob;
% second moment 
sen2m = zeros(1,nd);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);     end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);     end 
                end 
            else 
                % E[ON X ON X ON]           
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*triON(i,j,k);
            end 
        end 
    end 
 end  
sen2m = sen2m*jacob;
% variance 
varY2 = sum(cy(2:end).^2);
% mean
meanY2 = cy(1);
% constraint value 
c = 3*sqrt(varY2)-meanY2;
ceq = [];
% design sensitivities 
DC = ((3/2)*(1/sqrt(varY2))*(sen2m - 2*meanY2*sen1m)-sen1m)';
DCeq = [];

disp(c)
end % End of function 









            

















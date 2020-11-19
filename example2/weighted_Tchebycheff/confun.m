%% ========================================================================
% Example 1.1: Function of inequality constraint function and its grandient  
% (single-step approach)
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Two input required: 
% 1. design variables (dv) 
% 2. option for Gram matrix construction: 'QR' quadrature rule/ 'MC' Monte carlo integration (Quasi MC used)
%% ========================================================================
function [c, ceq, DC,  DCeq] = confun(dv)
global cntCon
%% Initialization
cntCon = cntCon + 1; % count the function call 
N = 2; % number of variables 
m = 1; % degree of ON (Orthonormal basis) for function y1
ms = 1; % degree of ON for score function
nd = 2;  % number of design variables 
% L_{N,m}
nA = nchoosek(N+m, m); % number of coefficients for function y0
nAs = nchoosek(N+ms, ms);% number of coefficients for score function 
A = max(nA, nAs);
% file name for loading data during 1st iteration 
FilNam = sprintf('gram.mat');
FilNam1 = sprintf('data.mat');
FilNam2 = sprintf('data2.mat');
% normalized mean (mu)
mu1 = 0; 
mu2 = 0;
mu = [mu1, mu2];
% coefficient of variation (sig)
sig1 = 0.4;
sig2 = 0.4;
sig = [sig1, sig2];
% correlation matrix  
rho12 = 0.4;
cor = [1,rho12; rho12,1];
% normalized mean's covariance matrix (cov)
cov = zeros(N,N);
for i=1:N
    for j=1:N
        if (i==j)
            cov(i,i) = sig(i)^2;
        else 
            cov(i,j) = rho12*sig(i)*sig(j);
        end 
    end 
end 
cs = zeros(nAs,nd);
load(FilNam);
load(FilNam1);

nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 1000000;
grSample = [nSampley, nSamples, nSampleo];
nSample = max(grSample); 
if (cntCon == 1)
    
% Update y1 and score function 
for L = 1:nSampley
        tmpY(L,1) = responY2(x(L,:),dv);
end 

infoMy = infoM([1:nSampley],[1:nA]);
infoMs = infoM([1:nSamples],[1:nAs]);
cy = (infoMy'*infoMy)\(infoMy'*tmpY);

%% Sensitivity analysis 
% first moment
sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs(j,i);
    end
end 
% scond moment 
sen2m = zeros(1,nd);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) ) 
                    sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);     end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);     end 
                end 
            else 
                % E[PsiXPsiXPsi]          
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*triON(i,j,k);
            end 
        end 
    end 
end  

varY2 = sum(cy(2:end).^2);
meanY2 = cy(1);

c = 3*sqrt(varY2)-meanY2;
ceq = [];

DC = ((3/2)*(1/sqrt(varY2))*(sen2m - 2*meanY2*sen1m)-sen1m)';
DCeq = [];
cy1 = cy;
disp(c)
save(FilNam2, 'cy1');

else 
    load(FilNam);
    load(FilNam1);
    load(FilNam2);
 for i = 1:nd
    x(:,i) = x(:,i) + dv(i)-dv0(i);
 end
x(nSample+1:end,:) = [];
 M = zeros(nSample,A); %monomial bases 
for i=1:A
    chkID = ID(i,:); %chkID is the same as how size of order  ex) X^2, X^3  
    nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
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

% Update y1 and score function 
for L = 1:nSampley 
        tmpY(L,1) = cy1'*QQ(1:nA,1:nA)*M(L,:)';
end 

infoMy = infoM([1:nSampley],[1:nA]);
infoMs = infoM([1:nSamples],[1:nAs]);
cy = (infoMy'*infoMy)\(infoMy'*tmpY);

%% Sensitivity analysis 
% first moment
sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs(j,i);
    end
end 
% scond moment 
sen2m = zeros(1,nd);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) ) 
                    sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);     end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);     end 
                end 
            else 
                % E[PsiXPsiXPsi]          
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*triON(i,j,k);
            end 
        end 
    end 
end  

varY2 = sum(cy(2:end).^2);
meanY2 = cy(1);

c = 3*sqrt(varY2)-meanY2;
ceq = [];

DC = ((3/2)*(1/sqrt(varY2))*(sen2m - 2*meanY2*sen1m)-sen1m)';
DCeq = [];

disp(c)


end % End of function 









            
















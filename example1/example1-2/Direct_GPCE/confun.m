%% ========================================================================
% Example 1.2: Function of inequality constraint function and its grandient  
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Two input required: 
% 1. design variables (dv) 
% 2. option for Gram matrix construction: 'QR' quadrature rule/ 'MC' Monte carlo integration (Quasi MC used)
%% ========================================================================
function [c, ceq, DC,  DCeq] = confun(dv)
global cntCon
double precision;
%% Initialization
cntCon = cntCon + 1; % count the function call 
N = 2; % number of variables 
m = 1; % degree of ON (Orthonormal basis) for function y1
ms = 1; % degree of ON for score function
nd = 2;  % number of design variables 
% L_{N,m}
nA = nchoosek(N+m, m); % number of coefficients for function y0
nAs = nchoosek(N+ms, ms);% number of coefficients for score function 
% file name for loading data during 1st iteration 
FilNam = sprintf('data.mat');
% normalized mean (mu)
mu1 = 1; 
mu2 = 1;
mu = [mu1, mu2];
% coefficient of variation (sig)
sig1 = 0.15;
sig2 = 0.15;
sig = [sig1, sig2];
% correlation matrix  
rho12 = -0.5;
cor = zeros(N,N);
% normalized mean's covariance matrix (cov)
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

%% load 'infoM', 'QQ', 'x', 'triON', 'ID', 'tmpS'
load(FilNam);
nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 1000000;
grSample = [nSampley, nSamples, nSampleo];
nSample = max(grSample); 
% Update y1 and score function 
for L = 1:nSample 
    if (L < nSampley +1) 
        tmpY(L,1) = responY2(x(L,:),dv);
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
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs1(j,i);
    end
end 
sen1m = sen1m*jacob;
% scond moment 
sen2m = zeros(1,nd);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs1(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs1(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs1(k,:);     end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs1(k,:);     end 
                end 
            else 
                % E[PsiXPsiXPsi]          
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs1(k,:)*triON(i,j,k);
            end 
        end 
    end 
 end  
sen2m = sen2m*jacob;

varY2 = sum(cy(2:end).^2);
meanY2 = cy(1);

c = 3*sqrt(varY2)-meanY2;
ceq = [];

DC = ((3/2)*(1/sqrt(varY2))*(sen2m - 2*meanY2*sen1m)-sen1m)';
DCeq = [];

% disp(c)
end % End of function 









            

















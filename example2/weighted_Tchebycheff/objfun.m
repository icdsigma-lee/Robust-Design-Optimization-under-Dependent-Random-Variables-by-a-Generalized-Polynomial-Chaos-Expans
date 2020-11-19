%% ========================================================================
% Example 1.1: Function of objective function and its grandient  (Single-step approach)
% Two input required: 
% 1. design variables (dv) 
% 2. option for Gram matrix construction: 'QR' quadrature rule/ 'MC' Monte carlo integration (Quasi MC used)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu)
%% ========================================================================
function [objValue, objGrad] = objfun(dv) 
global cntObj stat0 statf sopt w1 w2 
double precision;
%% Initialization
cntObj = cntObj + 1; % count the function call
N = 2; % number of variables  
m=4; % degree of ON (Orthonormal basis) for function y0   
ms = 1; % degree of ON for score function
nd = 2; % number of design variables 
% L_{N,m}
nA = nchoosek(N+m, m); % number of coefficients for function y0
nAs = nchoosek(N+ms, ms); % number of coefficients for score function 
grA = [nA, nAs];
A = max(grA);
% file name for saving data during 1st iteration 
FilNam = sprintf('gram.mat');
FilNam1 = sprintf('data.mat');
% zero mean (mu)
mu1 = 0; 
mu2 = 0;
mu = [mu1, mu2];
% coefficient of variation (sig)
sig1 = 0.4;
sig2 = 0.4;
sig = [sig1, sig2];
% correlation matrix  
rho12 = 0.4;
cor1 = zeros(N,N);
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
    cor1(i,j) = cov(i,j)/(sig(i)*sig(j));    
    end
end 

%% First iteration 

if (cntObj == 1) % generate Gram (G) and information matrix (infoM)
load(FilNam);


% Integration point number 
nGauss = ceil((m+1)/2)+10;
% Gauss points and weight values (Gaussian)
[xx, ww] =GaussHermite_2(nGauss);
xx1 = sqrt(2)*sig(1).*xx + mu(1);
xx2 = sqrt(2)*sig(2).*xx + mu(2);
ww1 = sqrt(1/pi).*ww;
ww2 = sqrt(1/pi).*ww;
tx = [xx1, xx2];
tw = [ww1, ww2];

% Generate normalized mean valued Gram-matrix 
%Cautions: max order: m=2
count = 0;


% Set sample size for 
% nSample for generic function 
% nSamples for score function 
% nSampleo for E[ON*ON*ON]

nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 1000000;
grSample = [nSampley, nSamples, nSampleo];
nSample = max(grSample); 
%% Expansion coefficient (Y)
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z1 = qrand(q,nSample);
x1 = norminv(z1,0,1);
ts = chol(cov,'lower'); % cov is normalized mean valued covariance matrix 
x = (ts*x1')';

for i=1:N
   x(:,i) = x(:,i) + mu(i);
end 
    
% least squares regression 
% information matrix
tmpY = zeros(nSampley,1);
tmpS = zeros(nSamples,nd);
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
onP = QQ*M';    
infoM = onP'; % information matrix: L X L_{N,m}

% (Part needed for updating at next steps) 
for L = 1:nSample 
    if (L < nSampley +1) 
        tmpY(L,1) = responY1(x(L,:),dv);
    end 
    if  (L <nSamples + 1) 
        scoredk = scoreGauss(x(L,:), mu, cov)';
        for p =1:nd
        tmpS(L,p) = scoredk(p);
        end  
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
% First moment of y0
sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs(j,i);
    end
end

% Second moment 
sen2m = zeros(1,nd);
triON = zeros(nA,nA,nAs);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
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
                % E[PsiXPsiXPsi] monte carlo integration 
                tmpYo = sum(infoM(:,i).*infoM(:,j).*infoM(:,k))/nSample;
                triON(i,j,k) = tmpYo;
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*tmpYo;
            end 
        end 
    end 
end  

varY1 = sum(cy(2:end).^2);
meanY1 = cy(1);

cy0 = cy;
dv0 = dv;
optMean = 4.4307; 
optStd = 1.1592;

scm = 31.5568;
scv = 17.0268;

crt1 = (w1/scm)*(meanY1 - optMean);
crt2 = (w2/scv)*(sqrt(varY1) - optStd);
objValue = [crt1,crt2];
objGrad = zeros(nd,2);
objGrad(:,1) = (w1/scm)*sen1m';
objGrad(:,2) = (w2/scv)*((1/2)*(1/sqrt(varY1))*(sen2m-2*meanY1*sen1m))'; 
save(FilNam1, 'infoM','x', 'triON','tmpS','cy0','cs','dv0')

else %(cntObj ~= 1)
    load(FilNam);
    load(FilNam1);
    nSampley = nA*3;
    nSamples = nAs*10;
    nSampleo = 1000000;
    grSample = [nSampley, nSamples, nSampleo];
    nSample = max(grSample); 
    % (Part needed for updating at next steps) 
    zTran = zeros(nd,nd);
for i = 1:nd
    x(:,i) = x(:,i) + dv(i)-dv0(i);
end 

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
 
for L = 1:nSample 
    if (L < nSampley +1) 
        tmpY(L,1) = cy0'*QQ*M(L,:)';
    end 
end 
infoMy = infoM([1:nSampley],[1:nA]);
infoMs = infoM([1:nSamples],[1:nAs]);
cy = (infoMy'*infoMy)\(infoMy'*tmpY);
%% Sensitivity analysis 
% First moment
sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs(j,i);
    end
end 

% Second moment 
sen2m = zeros(1,nd);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
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
                % E[PsiXPsiXPsi]
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*triON(i,j,k);
            end 
        end 
    end 
end 
 
varY1 = sum(cy(2:end).^2);
meanY1 = cy(1);

optMean = 4.4307; 
optStd = 1.1592;

scm = 31.5568;
scv = 17.0268;

crt1 = (w1/scm)*(meanY1 - optMean);
crt2 = (w2/scv)*(sqrt(varY1) - optStd);
objValue = [crt1,crt2];
objGrad = zeros(nd,2);
objGrad(:,1) = (w1/scm)*sen1m';
objGrad(:,2) = (w2/scv)*((1/2)*(1/sqrt(varY1))*(sen2m-2*meanY1*sen1m))'; 

end % End of function 

% record stat. info. at the initial and final design. 
switch sopt
    case 'pre'
        stat0 = [meanY1, sqrt(varY1)];
    case 'post'
        statf = [meanY1, sqrt(varY1)];
    otherwise 
end 




            

















%% ========================================================================
% Example 2: Function of objective function and its grandient  
% Input: design variables (dv) 
% Output: objective value and its design sensitivities
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [objValue, objGrad] = objfun(dv) 
global cntObj stat0 statf sopt w1 w2 
double precision;
%% Initialization
cntObj = cntObj + 1; % count function call #
N = 2; % # of random variables 
m=4; % order of GPCE for generic function 
ms = 1; % order of GPCE for score function
nd = 2; % # of design variables   
nA = nchoosek(N+m, m);% # of GPCE expansion coefficients for generic function
nAs = nchoosek(N+ms, ms); % # of GPCE expansion coefficients for score function 
grA = [nA, nAs];
A = max(grA);
% files for data at 1st iteration  
FilNam = sprintf('gram.mat');
FilNam1 = sprintf('data.mat');
% means vector (mu) 
mu1 = 0; 
mu2 = 0;
mu = [mu1, mu2];
% standard deviation (sig)
sig1 = 0.4;
sig2 = 0.4;
sig = [sig1, sig2];
% correlation matrix  
rho12 = 0.4;
cor1 = zeros(N,N);
cor = [1,rho12; rho12,1];
% covariance matrix (cov)
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
if (cntObj == 1) 
load(FilNam);
% Gauss point number 
nGauss = ceil((m+1)/2)+10;
% Gauss points and weight values (Gaussian)
[xx, ww] =GaussHermite_2(nGauss);
xx1 = sqrt(2)*sig(1).*xx + mu(1);
xx2 = sqrt(2)*sig(2).*xx + mu(2);
ww1 = sqrt(1/pi).*ww;
ww2 = sqrt(1/pi).*ww;
tx = [xx1, xx2];
tw = [ww1, ww2];

count = 0;

% nSampley for generic function 
% nSamples for score function 
% nSampleo for E[ON*ON*ON]

nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 1000000;
grSample = [nSampley, nSamples, nSampleo];
nSample = max(grSample); 

rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z1 = qrand(q,nSample);
x1 = norminv(z1,0,1);
ts = chol(cov,'lower'); 
x = (ts*x1')';

for i=1:N
   x(:,i) = x(:,i) + mu(i);
end 
    

tmpY = zeros(nSampley,1);
tmpS = zeros(nSamples,nd);
% Monomial basis (M) 
M = zeros(nSample,A); 
for i=1:A
    chkID = ID(i,:); %chkID: ex) (x1^(3),x2^(4))->(3,4) 
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
onP = QQ*M'; % whitening transformation (QQ=Wm in the paper)  
infoM = onP'; 

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
% First moment 
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
                % E[ON * ON * ON]
                tmpYo = sum(infoM(:,i).*infoM(:,j).*infoM(:,k))/nSample;
                triON(i,j,k) = tmpYo;
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*tmpYo;
            end 
        end 
    end 
end  
% variance
varY1 = sum(cy(2:end).^2);
% mean 
meanY1 = cy(1);

cy0 = cy;
dv0 = dv;


scm = 31.5568;
scv = 17.0268;

% objective value 
objValue =  w1*meanY1/scm + w2*sqrt(varY1)/scv;
% design sensitivities 
objGrad  = w1*(1/scm)*sen1m' + w2*((1/(2*scv))*(1/sqrt(varY1))*(sen2m - 2*meanY1*sen1m))';  

save(FilNam1, 'infoM','x', 'triON','tmpS','cy0','cs','dv0')

else %(cntObj ~= 1): the rest of the first iteration 
    load(FilNam);
    load(FilNam1);
    nSampley = nA*3;
    nSamples = nAs*10;
    nSampleo = 1000000;
    grSample = [nSampley, nSamples, nSampleo];
    nSample = max(grSample); 
    zTran = zeros(nd,nd);
for i = 1:nd
    x(:,i) = x(:,i) + dv(i)-dv0(i);
end 
% Monomial basis (M)
 M = zeros(nSample,A); 
for i=1:A
    chkID = ID(i,:); %chkID: ex) (x1^(2), x2^(3))-> (2,3)
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
                 % E[ON * ON * ON]
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*triON(i,j,k);
            end 
        end 
    end 
end 
% variance
varY1 = sum(cy(2:end).^2);
% mean
meanY1 = cy(1);


scm = 31.5568;
scv = 17.0268;

% objective value
objValue =  w1*meanY1/scm + w2*sqrt(varY1)/scv;
v
objGrad  = w1*(1/scm)*sen1m' + w2*((1/(2*scv))*(1/sqrt(varY1))*(sen2m - 2*meanY1*sen1m))';  

end % End of function 

% record state info. at initial and optimal design. 
switch sopt
    case 'pre'
        stat0 = [meanY1, sqrt(varY1)];
    case 'post'
        statf = [meanY1, sqrt(varY1)];
    otherwise 
end 




            

















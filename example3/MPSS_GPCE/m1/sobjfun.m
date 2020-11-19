%% ========================================================================
% Example 2: Function of objective function and its grandient (Single-Step GPCE)   
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Input required: design variables (dv) 
%% ========================================================================
function [objValue, objGrad] = sobjfun(dv)

%dv = ones(1,10)*30;
global cntObj stat0 statf sopt diffobj initobj
%ii = 0;
% tic
double precision;
w1 = 56.5744;
w2 = 17.0059;

%% Initialization
% number of variables
initobj = initobj + 1;
cntObj = cntObj + 1;
N = 7; 
m = 1; % ON degree for generic function 
ms = 1; % ON degree for score function
nd = 4; % design parameter size
N1 = 5;

% L_{N,m}
nA = nchoosek(N+m, m);
nAs = nchoosek(N+ms, ms);
A = max(nA,nAs);
nA1 = nchoosek(N1+m, m);
nAs1 = nchoosek(N1+ms, ms);

FilNam = sprintf('gram.mat');
FilNam1 = sprintf('data.mat');
% Define parameters of distributiosn 
mu5 = 10000; %kgm^3
mu6 = 2050000000; %Pa
mu7 = 200000; %N
% normalized mean vector for N variables 
mu = ones(1,N);
% transform to original mean vector 
muTr = [(dv(1)), (dv(2)), dv(3), dv(4), mu5, mu6, mu7];
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
%  The 4 by 4 covariance matrix for designing part

%% First iteration 
if (cntObj == 1)
load(FilNam); %% load ID, QQ 
if (initobj~=1)
    load(FilNam1);
end 
nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 2000000;
%nSampleo = nAo*nt;
grSample = [nSampley, nSamples, nSampleo];
nSample = max(grSample); 
% Sample generation for X1~X7
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z = qrand(q,nSample);
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

% Set sample size for 
% nSample for generic function 
% nSamples for score function 
% nSampleo for E[ON*ON*ON]

%% Expansion coefficient (Y)
% y1 function output (nSampleo by 1)
% standard least squares regression 
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
        M(:,i) = (x(:,nZeroID).^chkID(nZeroID));
    end 
    if (nZero == 2)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        M(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2));
    end 
    if (nZero == 3)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        M(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3));
    end 
    if (nZero == 4)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        M(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4));
    end 
end      

M0 = M;
M0(:,[INDEX0]) = [];
infoM0 = (QQ1*M0')';    

for L = 1:nSample 
    if (L < nSampley +1) 
        tmpY(L,1) = responY0(x(L,:),muTr);
    end 
    if  (L <nSamples + 1) 
        scoredk = scoreGauss(x(L,:), mu, cov,'full')';
        for p =1:nd
            tmpS(L,p) = scoredk(p);
        end  
    end 
end 
infoMy = infoM0([1:nSampley],:);
infoMs = infoM0([1:nSamples],[1:nAs1]);
infoMo = infoM0([1:nSampleo],:);

cs = zeros(nAs1,nd);
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
    for j = 1:min(nA1,nAs1)
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs(j,i);
    end
end 

sen1m = sen1m*jacob;
% second moment 
sen2m = zeros(1,nd);
if (initobj == 1)
triON = zeros(nA, nA, nAs);
end  
for i=1:nA1 %nA=m
        for j=1:nA1 %nA=m
            for k=1:nAs1 %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                end 
            else 
                   % E[PsiXPsiXPsi]
                index = [i,j,k];
				index1 = sort(index, 'descend');
				if (triON(index1(1), index1(2), index1(3)) == 0)
				triON(index1(1), index1(2), index1(3)) = sum(infoMo(:,i).*infoMo(:,j).*infoMo(:,k))/nSampleo;
                else 
				end 
                %tmpC = (infoMo'*infoMo)\(infoMo'*tmpYo);
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*triON(index1(1), index1(2), index1(3));
            end 
        end 
    end 
end 

sen2m = sen2m*jacob;

varY0 = sum(cy(2:end).^2);
meanY0 = cy(1);

objValue = 0.5*meanY0/w1 + 0.5*sqrt(varY0)/w2;
objGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varY0))*(sen2m - 2*meanY0*sen1m))'; 
disp(dv) 
disp(varY0)
disp(meanY0)

cy0 = cy;
dv0 = dv;

save(FilNam1, 'infoM0', 'x', 'triON', 'cs', 'M', 'cy0', 'dv0', '-v7.3');
% 2nd above iteration 
else %(cntObj~=1)
    load(FilNam);
    load(FilNam1);
    nSampley = nA*3;
    nSamples = nAs*10;
    nSampleo = 2000000;
    grSample = [nSampley, nSamples, nSampleo];
    nSample = max(grSample); 

    % (Part needed for updating at next steps) 
tmpY = zeros(nSampley,1);
zTran = zeros(nd, nd); 
for i = 1:nd
    zTran(i,i) = dv(i)/dv0(i);
end 
x(:,1:nd)= x(:,1:nd)*zTran; 
M1 = zeros(nSample,A);

for i=1:A
    chkID = ID(i,:); %chkID is the same as how size of order  ex) X^2, X^3  
    nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
    nZero = length(nZeroID);
    if (nZero == 0)
        M1(:,i) = 1;
    end 
    if (nZero == 1)
        M1(:,i) = (x(:,nZeroID).^chkID(nZeroID));
    end 
    if (nZero == 2)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        M1(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2));
    end 
    if (nZero == 3)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        M1(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3));
    end 
    if (nZero == 4)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        M1(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4));
    end 
end 

M1(:,[INDEX0]) = [];
tmpY = zeros(nSampley,1);
for L = 1:nSample 
    if (L < nSampley +1) 
        tmpY(L,1) = cy0'*QQ1*M1(L,:)';
    end 
end 

infoMy = infoM0([1:nSampley],:);
infoMs = infoM0([1:nSamples],[1:nAs1]);
infoMo = infoM0([1:nSampleo],:);

cy = (infoMy'*infoMy)\(infoMy'*tmpY);

 %% Sensitivity analysis 
% first moment
 jacob = zeros(nd,nd);
for i = 1:nd
    jacob(i,i) = 1/dv(i);
end 

sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA1,nAs1)
    sen1m(1,i) = sen1m(1,i) + cy(j)*cs(j,i);
    end
end 

sen1m = sen1m*jacob;

% second moment 
sen2m = zeros(1,nd);
    for i=1:nA1 %nA=m
        for j=1:nA1 %nA=m
            for k=1:nAs1 %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  cy(i)*cy(j)*cs(k,:);  end 
                end 
            else 
                % E[PsiXPsiXPsi]
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs(k,:)*triON(i,j,k);
            end 
        end 
    end 
    end     

sen2m = sen2m*jacob;
    
varY0 = sum(cy(2:end).^2);
meanY0 = cy(1);

objValue = 0.5*meanY0/w1 + 0.5*sqrt(varY0)/w2;
objGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varY0))*(sen2m - 2*meanY0*sen1m))'; 
disp('dv:')
disp(dv) 
disp('varY0:')
disp(varY0)
disp('meanY0:')
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



            

















%% ========================================================================
% Example 1.2: Function of objective function and its grandient  
% Input: design variables (dv) 
% Output: objective value and its design sensitivities
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [objValue, objGrad] = objfun(dv) 
global cntObj 
double precision;
%% Initialization
cntObj = cntObj + 1; % count function call #
opt = 'QR';
N = 2; % # of random variables 
m=4; % order of GPCE for generic function 
ms = 1; % order of GPCE for score function
nd = 2; % # of design variables  
nA = nchoosek(N+m, m); % # of GPCE expansion coefficients for generic function
nAs = nchoosek(N+ms, ms); % # of GPCE expansion coefficients for score function 
% files for data at 1st iteration 
FilNam = sprintf('data.mat');
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

%% First iteration 
if (cntObj == 1) 

% Graded lexicographical ordered indeterminates (ID)
cnt = 0;
for m0=1:max(m,ms)+1
    mm = m0-1; %total degree 
for i1=mm+1:-1:1
    for i2=mm+1:-1:1               
                             j1 = i1-1;
                             j2 = i2-1;
                             if ((j1+j2)==mm)
                                 cnt = cnt + 1;
                                 ID(cnt,:) = [j1 j2];
                             end
                        
    end 
end 
end              
disp('Compeletion of ID(graded lexicographical order)')

switch opt
    case 'QR' %Option for quadrature 
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

% Monomial moment matrix (G) 
grA = [nA, nAs];
A = max(grA);
G = zeros(A,A); 
idm = 2*m+1; 
Mom = zeros(idm, idm); % moments-data matrix 
for iRow=1:A 
    for iCol=iRow:A 
        chkID = ID(iRow,:) + ID(iCol,:); 
        nZeroID = find(chkID~=0);
        nZero = length(nZeroID);  
        tmp = 0; 
        tmp2 = Mom(chkID(1)+1, chkID(2)+1);
       if (tmp2 == 0)
        if (nZero == 0)
            tmp = 1;
        end 
        if (nZero == 1) 
            for k = 1:nGauss 
                tmp = tmp + bas(tx(k,nZeroID),chkID(nZeroID))*tw(k,nZeroID);
            end 
        end 
        if (nZero == 2)
            id1 = nZeroID(1); id2 = nZeroID(2);
            tmpMu2 = [mu(id1); mu(id2)];
            cov2 = [sig(id1)^2, cor(id1,id2)*sig(id1)*sig(id2); cor(id1,id2)*sig(id1)*sig(id2), sig(id2)^2];         
            for k1 = 1:nGauss
                for k2 = 1:nGauss  
                    biX = [tx(k1,id1);tx(k2,id2)];                                                     
                    biTrans = nGaussDen(biX, tmpMu2, cov2)/(nGaussDen(tx(k1,id1), mu(id1), sig(id1)^2)*nGaussDen(tx(k2,id2), mu(id2), sig(id2)^2));
                    tmp = tmp + bas(tx(k1,id1),chkID(id1))*bas(tx(k2,id2),chkID(id2))*biTrans*tw(k1,id1)*tw(k2,id2);                    
                end 
            end 
        end 
        G(iRow, iCol) = tmp;
        Mom(chkID(1)+1, chkID(2)+1) = tmp;
       else 
           G(iRow, iCol) = tmp2;
       end 
    end 
end 

    case 'MC' % Option for sampling method  
nSample = 500000;
% Quasi-Monte Carlo sampling 
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen'); 
q = qrandstream(p);
z = qrand(q,nSample);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);
T12 = chol(cov,'lower');
x(:,1:2) = (T12*x(:,1:2)')';
for i=1:2
    x(:,i) = x(:,i) + mu(i);
end 
count = 0;
grA = [nA, nAs];
A = max(grA);
G = zeros(A,A); 
for iRow=1:A 
    for iCol=iRow:A 
        chkID = ID(iRow,:) + ID(iCol,:); 
        nZeroID = find(chkID~=0); 
        nZero = length(nZeroID);  
        tmp = 0;
        if (nZero == 0)
            tmp = 1;
        end 
        if (nZero == 1) 
                tmp = sum(x(:,nZeroID).^chkID(nZeroID))/nSample;
        end 
        if (nZero == 2)
            id1 = nZeroID(1); id2 = nZeroID(2);                          
                tmp = sum ((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)))/nSample;                    
        end 
        count = count + 1;
        G(iRow, iCol) = tmp;
    end 
end

    otherwise 
    disp('wrong option was used for monomial moment matrix')
end % end of the construction of monomial moment matrix  
        
for iRow=1:A 
    for iCol=iRow+1:A 
        if (abs(G(iRow,iCol)) < 1E-15)
        G(iRow,iCol) = 0;
        end 
        G(iCol, iRow) = G(iRow,iCol);
    end 
end 
% Cholesky decomposition of monomial moment matrix 
Q = chol(G,'lower');
% whitening transformation matrix (QQ)  
QQ = inv(Q);

% nSample for generic function 
% nSamples for score function 
% nSampleo for E[ON*ON*ON]

nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 1000000;
grSample = [nSampley, nSamples, nSampleo];
nSample = max(grSample); 

if (opt == 'QR') 
% Note that herein x is presented as z in the theorem, that is transformed
% variables, and its mean value vector is ones in scaling transformation  
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
end 
x([nSample+1:end],:) = [];
tmpY = zeros(nSampley,1);
tmpS = zeros(nSamples,nd);
M = zeros(nSample,A); 
for i=1:A
    chkID = ID(i,:);   %chkID: ex) (x1^(1),x2^(2))->(1,2)
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
onP = QQ*M';    
infoM = onP'; % information matrix

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
cs1 = zeros(nAs,nd);
cy = (infoMy'*infoMy)\(infoMy'*tmpY);
for i=1:nd 
    cs1(:,i) = (infoMs'*infoMs)\(infoMs'*tmpS(:,i));
end 
cs1(1,:) = 0;
%% Sensitivity analysis 
% First moment 
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

% Second moment 
sen2m = zeros(1,nd);
triON = zeros(nA,nA,nAs);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1))
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
                % E[ON X ON X ON]
                tmpYo = sum(infoM(:,i).*infoM(:,j).*infoM(:,k))/nSample;
                triON(i,j,k) = tmpYo;
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs1(k,:)*tmpYo;
            end 
        end 
    end 
end  


sen2m = sen2m*jacob;
% variance
varY1 = sum(cy(2:end).^2);
% mean
meanY1 = cy(1);
% objective value 
objValue = sqrt(varY1)/45;
% design sensitivities
objGrad = ((1/90)*(1/sqrt(varY1))*(sen2m-2*meanY1*sen1m))'; 
save(FilNam, 'infoM', 'QQ', 'x', 'triON', 'ID', 'cs1')

else %(cntObj ~= 1): the rest of the first iteration 
    load(FilNam); 
    nSampley = nA*3;
    nSamples = nAs*10;
    nSampleo = 1000000;
    grSample = [nSampley, nSamples, nSampleo];
    nSample = max(grSample); 
for L = 1:nSample 
    if (L < nSampley +1) 
        tmpY(L,1) = responY1(x(L,:),dv);
    end 
end 
infoMy = infoM([1:nSampley],[1:nA]);
infoMs = infoM([1:nSamples],[1:nAs]);

cy = (infoMy'*infoMy)\(infoMy'*tmpY);

%% Sensitivity analysis 
% First moment
 jacob = zeros(nd,nd);
for i = 1:nd
    jacob(i,i) = 1/dv(i);
end 

sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) +cy(j)*cs1(j,i);
    end
end 

sen1m = sen1m*jacob;

% Second moment 
sen2m = zeros(1,nd);
for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) 
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
                %  E[ON X ON X ON]
                sen2m(1,:) = sen2m(1,:) + cy(i)*cy(j)*cs1(k,:)*triON(i,j,k);
            end 
        end 
    end 
end 

sen2m = sen2m*jacob;
% variance 
varY1 = sum(cy(2:end).^2);
% mean
meanY1 = cy(1);

objValue = sqrt(varY1)/45;
objGrad = ((1/90)*(1/sqrt(varY1))*(sen2m-2*meanY1*sen1m))'; 
display(meanY1);
display(varY1);
end % End of function 









            

















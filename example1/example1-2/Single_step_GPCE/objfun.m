%% ========================================================================
% Example 1.2: Function of objective function and its grandient  (Single-step approach)
% Two input required: 
% 1. design variables (dv) 
% 2. option for Gram matrix construction: 'QR' quadrature rule/ 'MC' Monte carlo integration (Quasi MC used)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu)
%% ========================================================================
function [objValue, objGrad] = objfun(dv) 
global cntObj 
double precision;
%% Initialization
cntObj = cntObj + 1; % count the function call
opt = 'QR';
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
% save 'Gram matrix', 'expectation of triple products of orthogonal
% polynomials 
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
%% First iteration 

if (cntObj == 1) % generate an monmial moment matrix (G) 

% Generate a set of graded lexicographical ordered indeterminates (ID)
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
    case 'QR' %quadrature 
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

% generate ones means-valued monomial moment-matrix 

G = zeros(A,A); % Initialize gram matrix 
idm = 2*m+1; % max. degree of monomial moments in a two-dimensional monomial moment matrix  
Mom = zeros(idm, idm); % moments-data matrix 
for iRow=1:A 
    for iCol=iRow:A 
        % Estimate index of E(X_row*X_col) 
        chkID = ID(iRow,:) + ID(iCol,:); %chkID is the same as how size of order  ex) X^2, X^3  
        nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
        nZero = length(nZeroID);  
        tmp = 0; %initialization (Important)
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
            % Set probabilistic property
            id1 = nZeroID(1); id2 = nZeroID(2);
            tmpMu2 = [mu(id1); mu(id2)];
            % Generate covariance matrix 
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
        %disp(count);
        %disp('# of Gauss:'); disp(nGauss);
        %disp('# of variable:'); disp(nZero);
    end 
end 

    case 'MC' % Sampling method for Gram matrix  

nSample = 1000000;
% Quasi-Monte Carlo sampling 
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen'); 
q = qrandstream(p);
z = qrand(q,nSample);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);
% Transformation to correlated random variables (x) 
T12 = chol(cov,'lower');
x(:,1:2) = (T12*x(:,1:2)')';
for i=1:2
    x(:,i) = x(:,i) + mu(i);
end 
%% Generate normalized mean valued Gram-matrix 
%Cautions: max order: m=2
count = 0;
grA = [nA, nAs];
A = max(grA);
G = zeros(A,A); % Initialize gram matrix 
for iRow=1:A 
    for iCol=iRow:A 
        % Estimate index of E(X_row*X_col) 
        chkID = ID(iRow,:) + ID(iCol,:); %chkID is the same as how size of order  ex) X^2, X^3  
        nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
        nZero = length(nZeroID);  
        tmp = 0; %initialization (Important)
        if (nZero == 0)
            tmp = 1;
        end 
        if (nZero == 1) 
                tmp = sum(x(:,nZeroID).^chkID(nZeroID))/nSample;
        end 
        if (nZero == 2)
            % Set probabilistic property
            id1 = nZeroID(1); id2 = nZeroID(2);                          
                tmp = sum ((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)))/nSample;                    
        end 
        count = count + 1;
        G(iRow, iCol) = tmp;
        %disp(count);
        %disp('# of Gauss:'); disp(nGauss);
        %disp('# of variable:'); disp(nZero);
    end 
end

    otherwise 
    disp('wrong option was wrong for Gram matrix generation')
end % end of Gram matrix generation using selected method 
        
for iRow=1:A 
    for iCol=iRow+1:A 
        if (abs(G(iRow,iCol)) < 1E-15)
        G(iRow,iCol) = 0;
        end 
        G(iCol, iRow) = G(iRow,iCol);
    end 
end 
% Cholesky decomposition of Gram matrix
Q = chol(G,'lower');
% ON transformation matrix 
QQ = inv(Q);

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
if (opt == 'QR') % if opt == MC, samples are reused. 
% Note that herein x is presented as z in the theorem, that is transformed
% variable, so its mean value vector is ones 
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
end 
x([nSample+1:end],:) = [];
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

sen2m = sen2m*jacob;

varY1 = sum(cy(2:end).^2);
meanY1 = cy(1);
objValue = sqrt(varY1)/45;
cy0 = cy;
dv0 = dv;
objGrad = ((1/90)*(1/sqrt(varY1))*(sen2m-2*meanY1*sen1m))'; 
save(FilNam, 'infoM', 'QQ', 'x', 'triON', 'ID', 'tmpS', 'cy0', 'cs', 'dv0')

else %(cntObj ~= 1)
    load(FilNam); 
    nSampley = nA*3;
    nSamples = nAs*10;
    nSampleo = 1000000;
    grSample = [nSampley, nSamples, nSampleo];
    nSample = max(grSample); 
    % (Part needed for updating at next steps)
    % transformation of z to z'
    zTran = zeros(nd,nd);
for i = 1:nd
    zTran(i,i) = dv(i)/dv0(i);
end 
 x = x*zTran;    
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
    jacob = zeros(nd,nd);
for i = 1:nd
    jacob(i,i) = 1/dv(i);
end 
sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) +cy(j)*cs(j,i);
    end
end 

sen1m = sen1m*jacob;

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

sen2m = sen2m*jacob;

varY1 = sum(cy(2:end).^2);
meanY1 = cy(1);

objValue = sqrt(varY1)/45;
objGrad = ((1/90)*(1/sqrt(varY1))*(sen2m-2*meanY1*sen1m))'; 
display(meanY1);
display(varY1);
end % End of function 









            

















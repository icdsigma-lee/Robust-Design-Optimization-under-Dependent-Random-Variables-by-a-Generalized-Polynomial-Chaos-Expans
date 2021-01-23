%% ========================================================================
%  function to generate monomial moment matrix 
%  output: gram.mat
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function genGramMatrix  
N = 2; % # of random variables 
m=4; % order of GPCE for generic function 
ms = 1; % order of GPCE for score function
nd = 2; % # of design variables  
nA = nchoosek(N+m, m); % # of GPCE expansion coefficients for generic function 
nAs = nchoosek(N+ms, ms);  % # of GPCE expansion coefficients for score function 
% files for data at 1st iteration 
FilNam = sprintf('gram.mat');
opt = 'QR'; % option for quadrature 
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
% original covariance matrix (covT)
covT = zeros(N,N);
for i=1:N
    for j=1:N
        if (i==j)
            covT(i,i) = (sig(i))^2;
        else 
            covT(i,j) = rho12*sig(i)*sig(j);
        end 
    end 
end

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

count = 0;
grA = [nA, nAs];
A = max(grA);
% monomial moment matrix (G=Gm in paper)
G = zeros(A,A); 
for iRow=1:A 
    for iCol=iRow:A  
        chkID = ID(iRow,:) + ID(iCol,:); %chkID: ex) (row)(x1^(1),x2^(2))*(column)(x1^(2),x2^(3))->(3,5)
        nZeroID = find(chkID~=0); 
        nZero = length(nZeroID);  
        tmp = 0; 
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
        count = count + 1;
        G(iRow, iCol) = tmp;
    end 
end 

    case 'MC' % Option for Sampling method 

nSample = 1000000;
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
% Monomial moment matrix (G)
G = zeros(A,A); 
for iRow=1:A 
    for iCol=iRow:A 
        chkID = ID(iRow,:) + ID(iCol,:); %chkID: ex) (row)(x1^(1),x2^(2))*(column)(x1^(2),x2^(3))->(3,5)
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
    disp('wrong option was selected for monomial moment matrix')
end % end of the construction of monomial moment matrix 
        
for iRow=1:A 
    for iCol=iRow+1:A 
        if (abs(G(iRow,iCol)) < 1E-15)
        %G(iRow,iCol) = 0;
        end 
        G(iCol, iRow) = G(iRow,iCol);
    end 
end 
% Cholesky decomposition 
Q = chol(G,'lower');
% whitening transformation matrix (QQ=Wm in the paper) 
QQ = inv(Q);

save(FilNam, 'ID', 'QQ');
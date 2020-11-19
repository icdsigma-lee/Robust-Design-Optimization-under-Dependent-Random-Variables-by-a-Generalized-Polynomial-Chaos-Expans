%% ========================================================================
%  function to generate monomial moment matrix or gram matrix (GPCE)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function genGramMatrix
% user defined initialization 
N = 2; % number of variables  
m=4; % degree of ON (Orthonormal basis) for function y0   
% L_{N,m}
opt = 'QR';
nA = nchoosek(N+m, m); % number of coefficients for function y0
% file name for saving data during 1st iteration 
FilNam = sprintf('gram.mat');
% zero mean (mu)
mu1 = 1; 
mu2 = 1;
mu = [mu1, mu2];
% coefficient of variation (sig)
sig1 = 0.15;
sig2 = 0.15;
sig = [sig1, sig2];
% correlation matrix  
rho12 = -0.5;
cor= zeros(N,N);
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

% Generate a set of graded lexicographical ordered indeterminates (ID)
cnt = 0;
for m0=1:m+1
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
% Generate normalized mean valued Gram-matrix 
%Cautions: max order: m=2

G = zeros(nA,nA); % Initialize gram matrix 
idm = 2*m+1; % max. degree of monomial moments in a two-dimensional monomial moment matrix  
Mom = zeros(idm, idm); % moments-data matrix  
for iRow=1:nA 
    for iCol=iRow:nA 
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
                tmp = tmp + (tx(k,nZeroID)^chkID(nZeroID))*tw(k,nZeroID);
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
                    tmp = tmp + (tx(k1,id1)^chkID(id1))*(tx(k2,id2)^chkID(id2))*biTrans*tw(k1,id1)*tw(k2,id2);                    
                end 
            end 
        end
            G(iRow, iCol) = tmp;
            Mom(chkID(1)+1, chkID(2)+1)=tmp;
     else
            G(iRow, iCol) = tmp2;
     end 
        %disp(count);
        %disp('# of Gauss:'); disp(nGauss);
        %disp('# of variable:'); disp(nZero);
    end 
end
    case 'MC'
        % on working 
    otherwise 
end 

for iRow=1:nA 
    for iCol=iRow+1:nA 
        if (abs(G(iRow,iCol)) < 1E-15)
        %G(iRow,iCol) = 0;
        end 
        G(iCol, iRow) = G(iRow,iCol);
    end 
end 
Q = chol(G,'lower');
QQ = inv(Q);
save(FilNam, 'ID', 'QQ');
%% ========================================================================
%  function to generate monomial moment matrix or gram matrix (GPCE)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% output : monomial moment matrix (N=7)
%% ========================================================================

function genGramMatrix

% user defined initialization 
N = 7; % number of variables  
m = 3; % ON degree for generic function 
ms = 1; % ON degree for score function
mo = 2; % ON dgree for E[ON*ON*ON]   
nd = 4; % design parameter size 
% L_{N,m}
nA = nchoosek(N+m, m);
nAs = nchoosek(N+ms, ms);
nAo = nchoosek(N+mo, mo);
opt = 'MC';
FilNam = sprintf('gram.mat');
% zero mean (mu)
mu = ones(1,N);
% Define parameters of distributiosn 
mu5 = 10000; %kgm^3
mu6 = 2050000000; %Pa
mu7 = 200000; %N
% normalized mean vector for N variables 
mu = ones(1,N);
% transform to original mean vector 
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


nSample = 5000000;
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

%% Generate a set of graded lexicographical ordered indeterminates (ID)
cnt = 0;
for m0=1:m+1
    mm = m0-1; %total degree 
for i1=mm+1:-1:1
    for i2=mm+1:-1:1              
        for i3=mm+1:-1:1    
            for i4=mm+1:-1:1    
                for i5=mm+1:-1:1    
                	for i6=mm+1:-1:1   
                		for i7=mm+1:-1:1   
                            j1 = i1-1;
                            j2 = i2-1;
                            j3 = i3-1;
                            j4 = i4-1;
                            j5 = i5-1;
                            j6 = i6-1;
                            j7 = i7-1;
            if ((j1+j2+j3+j4+j5+j6+j7)==mm)
                cnt = cnt + 1;
                ID(cnt,:) = [j1 j2 j3 j4 j5 j6 j7];
            end 
                        end
                    end
                end 
         	end            
        end            
    end 
end 
end       
disp('Compeletion of ID(graded lexicographical order)')
switch opt
    case 'QR' %quadrature 
    case 'MC'
%% Generate Gram-matrix 
%Cautions: max order: m=2
count = 0;
grA = [nA, nAs];
A = max(grA);
G = zeros(A,A); % Initialize gram matrix 
idm = 2*m + 1;
Mom = zeros(idm, idm, idm, idm, idm, idm, idm); % moments-data matrix 
for iRow=1:A 
    for iCol=iRow:A 
        % Estimate index of E(X_row*X_col) 
        chkID = ID(iRow,:) + ID(iCol,:); %chkID is the same as how size of order  ex) X^2, X^3  
        nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
        nZero = length(nZeroID);  
        tmp = 0; %initialization (Important)
        tmp2 = Mom(chkID(1)+1, chkID(2)+1, chkID(3)+1, chkID(4)+1, chkID(5)+1, chkID(6)+1, chkID(7)+1); 
        if (tmp2 == 0)
            if (nZero == 0)
                tmp = 1;
            end 
            if (nZero == 1) 
                    tmp =  sum(x(:,nZeroID).^chkID(nZeroID))/nSample;
            end 
            if (nZero == 2)
            % Set probabilistic property
                id1 = nZeroID(1); id2 = nZeroID(2);
            % Generate covariance matrix 
                 tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)))/nSample;
            end 
             if (nZero == 3)
            % Set probabilistic property
               id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3);
            % Generate covariance matrix 
               tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)))/nSample;
            end 
             if (nZero == 4)
               % Set probabilistic property
               id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4);
               % Generate covariance matrix 
               tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)))/nSample;
             end 
            if (nZero == 5)
             % Set probabilistic property
                id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);
                    % Generate covariance matrix 
               tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)))/nSample;
            end 
             if (nZero == 6)
                % Set probabilistic property
                id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6);
                % Generate covariance matrix 
                tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)))/nSample;
             end 
             if (nZero == 7)
                % Set probabilistic property
                 id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7);
                 % Generate covariance matrix 
                tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)))/nSample;
             end 
             count = count + 1;
             G(iRow, iCol) = tmp;
             Mom(chkID(1)+1, chkID(2)+1, chkID(3)+1, chkID(4)+1, chkID(5)+1, chkID(6)+1, chkID(7)+1) = tmp;
         %disp(count);
         %disp('# of Gauss:'); disp(nGauss);
            %disp('# of variable:'); disp(nZero);
        else 
            G(iRow, iCol) = tmp2;
        end 
     end 
end 

% tGram = toc;
% recover from indefinite matrix
% [V,D,W] = eig(G);
% for i=1:nA 
%     if ((D(i,i) <0) || (D(i,i) == 0))
%         D(i,i) = 1E-13;
%     else end
% end 
% G1 = (V*D*W');
% Cholesky decomposition of Gram matrix
% Set sample size for 
% nSample for generic function 
% nSamples for score function 
% nSampleo for E[ON*ON*ON]
        % on working 
    otherwise 
end 

for iRow=1:A 
    for iCol=iRow+1:A 
        if (abs(G(iRow,iCol)) < 1E-15)
        G(iRow,iCol) = 0;
        end 
        G(iCol, iRow) = G(iRow,iCol);
    end 
end 

[V,D,W] = eig(G);
for i=1:nA 
    if ((D(i,i) <0) || (D(i,i) == 0))
        D(i,i) = 1E-13;
    else end
end 
G1 = (V*D*W');

Q = chol(G1,'lower');
QQ = inv(Q);
save(FilNam, 'ID', 'QQ', 'G1');
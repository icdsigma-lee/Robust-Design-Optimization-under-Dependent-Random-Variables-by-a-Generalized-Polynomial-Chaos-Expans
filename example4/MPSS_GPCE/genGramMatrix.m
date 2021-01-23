%% ========================================================================
%  function to generate monomial moment matrix or gram matrix (GPCE)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% output : monomial moment matrix (N=10)
%% ========================================================================
function genGramMatrix
tic
% user defined initialization 
N = 10; % number of random variables  
m = 2; % ON degree for generic function 
ms = 3; % ON degree for score function
nd = 10; % design parameter size 
% L_{N,m}
nA = nchoosek(N+m, m);
nAs = nchoosek(N+ms, ms);
A = max(nA, nAs);
m = max(ms,m); 

opt = 'MC';
FilNam = sprintf('frtstm3cor2.mat'); % data saving for the first step 

% normalized mean vector for N variables 
NMN = ones(1,N);

% standard deviatiaon from
cov= 0.05; % coeff. of variation 
SIG = [cov, cov, cov, cov, cov, cov, cov, cov, cov, cov];

% log-normal  
for i = 1:nd
mpar(i) = log(NMN(i)^2/sqrt(SIG(i)^2 + NMN(i)^2)); % paramter mu 
spar(i) = sqrt(log(1 + ((SIG(i))^2)/(NMN(i)^2))); % paramter sig
end

% Correlation coefficient matrix in normal dist. 
cor = zeros(N,N);
for i=1:N
    for j=i:N
        if (i~=j)
            cor(i,j) = 0.2;
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
% Covariance matrix for X1~X10
cova = zeros(nd,nd); 
for i =1:nd
    for j = i:nd 
        cova(i,j) = cor(i,j)*spar(i)*spar(j);
    end 
end 
for i=1:nd 
    for j=i:nd 
        cova(j,i) = cova(i,j);
    end 
end

nSample = 5000000;
%% Sampling
% Sample generation for X1~X10
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z = qrand(q,nSample);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);

% transformation for x1~x8
trfo = chol(cova,'lower');
x(:,1:nd) = (trfo*x(:,1:nd)')';

% lognormal 
for i=1:nd
    x(:,i) = exp(x(:,i) + mpar(i));
end 

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
                            for i8=mm+1:-1:1
                                for i9=mm+1:-1:1
                                    for i10=mm+1:-1:1
                            j1 = i1-1;
                            j2 = i2-1;
                            j3 = i3-1;
                            j4 = i4-1;
                            j5 = i5-1;
                            j6 = i6-1;
                            j7 = i7-1;
                            j8 = i8-1;
                            j9 = i9-1;
                            j10 = i10-1;
            if ((j1+j2+j3+j4+j5+j6+j7+j8+j9+j10)==mm)
                cnt = cnt + 1;
                ID(cnt,:) = [j1 j2 j3 j4 j5 j6 j7 j8 j9 j10];
            end 
                                        end
                                    end
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
        % under construction 
    case 'MC'
%% Generate Gram-matrix 
%Cautions: max order: m=2
count = 0;
GRM = zeros(nA, nA); % Initialize GRAM matrix 
idm = 2*m + 1;
momtrix = zeros(idm, idm, idm, idm, idm, idm, idm, idm, idm, idm); % moments-data matrix 
for iRow=1:A 
    for iCol=iRow:A 
        % Estimate index of E(X_row*X_col) 
        chkID = ID(iRow,:) + ID(iCol,:); %chkID is the same as how size of order  ex) X^2, X^3  
        nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
        nZero = length(nZeroID);  
        tmp = 0; %initialization (Important)
        tmp2 = momtrix(chkID(1)+1, chkID(2)+1, chkID(3)+1, chkID(4)+1, chkID(5)+1, chkID(6)+1, chkID(7)+1, chkID(8)+1, chkID(9)+1, chkID(10)+1); 
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
             if (nZero == 8)
                % Set probabilistic property
                 id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7); id8 = nZeroID(8); 
                 % Generate covariance matrix 
                tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)).*(x(:,id8).^chkID(id8)))/nSample;
             end 
             if (nZero == 9)
                % Set probabilistic property
                 id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7); id8 = nZeroID(8); id9 = nZeroID(9); 
                 % Generate covariance matrix 
                tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)).*(x(:,id8).^chkID(id8)).*(x(:,id9).^chkID(id9)))/nSample;
              end 
             if (nZero == 10)
                % Set probabilistic property
                 id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7); id8 = nZeroID(8); id9 = nZeroID(9);  id10 = nZeroID(10); 
                 % Generate covariance matrix 
                tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)).*(x(:,id8).^chkID(id8)).*(x(:,id9).^chkID(id9)).*(x(:,id10).^chkID(id10)))/nSample;
              end 
             if (nZero == 11)
                % Set probabilistic property
                 id1 = nZeroID(1); id2 = nZeroID(2); id3 = nZeroID(3); id4 = nZeroID(4); id5 = nZeroID(5);  id6 = nZeroID(6); id7 = nZeroID(7); id8 = nZeroID(8); id9 = nZeroID(9);  id10 = nZeroID(10);   id11 = nZeroID(11); 
                 % Generate covariance matrix 
                tmp = sum((x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6)).*(x(:,id7).^chkID(id7)).*(x(:,id8).^chkID(id8)).*(x(:,id9).^chkID(id9)).*(x(:,id10).^chkID(id10)).*(x(:,id11).^chkID(id11)))/nSample;
              end 
             count = count + 1;
             GRM(iRow, iCol) = tmp;
             momtrix(chkID(1)+1, chkID(2)+1, chkID(3)+1, chkID(4)+1, chkID(5)+1, chkID(6)+1, chkID(7)+1, chkID(8)+1, chkID(9)+1, chkID(10)+1) = tmp;
        else 
            GRM(iRow, iCol) = tmp2;
        end 
     end 
end 

    otherwise 
end 

for iRow=1:A 
    for iCol=iRow+1:A 
        if (abs(GRM(iRow,iCol)) < 1E-15)
        GRM(iRow,iCol) = 0;
        end 
        GRM(iCol, iRow) = GRM(iRow,iCol);
    end 
end 
% Repair ill-conditioning of GRM 
[V,D,W] = eig(GRM);
for i=1:nA 
    if ((D(i,i) <0) || (D(i,i) == 0))
        D(i,i) = 1E-13;
    else end
end 
GRM1 = (V*D*W');

Q = chol(GRM1,'lower');
ORN = inv(Q); % orthonormalization matrix  
save(FilNam, 'ID','ORN');
toc
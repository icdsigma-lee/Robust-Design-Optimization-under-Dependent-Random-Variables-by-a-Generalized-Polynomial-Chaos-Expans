%% ========================================================================
% Example 2: Function of inequality constraint function and its grandient (Direct GPCE)   
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Input required: design variables (dv) 
%% ========================================================================
function[c, ceq, DC,  DCeq] = confun(dv)

global cntCon con
double precision;
%% Initialization
% number of variables
cntCon = cntCon + 1;
N = 7; 
m =3; % ON degree for generic function 
ms = 1; % ON degree for score function
nd = 4; % design parameter size 
% L_{N,m}
nA = nchoosek(N+m, m);
nAs = nchoosek(N+ms, ms);
% nA for y1
nA1 = nchoosek(N-2+m, m);
nAs1 = nchoosek(N-2+ms, ms);
% nA for y2
nA2 = nchoosek(N-2+m, m);
nAs2 = nchoosek(N-2+ms, ms);

A = max(nA, nAs);
FilNam = sprintf('gram.mat');
FilNam1 = sprintf('data.mat');
FilNam2 = sprintf('data2.mat');
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
cov = zeros(nd,nd);
cov(1:2,1:2) = cov12;
cov(3:4,3:4) = cov34;

% load 
load(FilNam);   
load(FilNam1);


% Set sample size for 
% nSample for generic function 
% nSamples for score function 
% nSampleo for E[ON*ON*ON]

nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 2000000;
grSample = [nSampley, nSamples, nSampleo];
nSample = max(grSample); 
if (cntCon == 1) 
% create triple product of O.G. 
%% Expansion coefficient (Y)
% y1 function output (nSampleo by 1)
tmpY1 = zeros(nSampley,1);
tmpY2 = zeros(nSampley,1);
for L=1:nSampley
        if (L < nSampley +1) 
            tmpY1(L,1) = responY1(x(L,:), muTr);
            tmpY2(L,1) = responY2(x(L,:), muTr);
        end 
end 

tmpS1 = zeros(nSamples, nd-1);
tmpS2 = zeros(nSamples, nd-1);

for L = 1:nSample 
    if  (L <nSamples + 1) 
        scoredk1 = scoreGauss(x(L,:), mu, cov,'sub1')';
        scoredk2 = scoreGauss(x(L,:), mu, cov,'sub2')';
        for p =1:3
            tmpS1(L,p) = scoredk1(p);
            tmpS2(L,p) = scoredk2(p);
        end  
    end 
end 
% monomial basis for y1 and y2
M1 = M;
M2 = M;
M1(:,[INDEX1]) = [];
M2(:,[INDEX2]) = [];
infoM1 = (QQ2*M1')';
infoM2 = (QQ3*M2')';


infoMo1 = infoM1([1:nSampleo],:);
infoMo2  = infoM2([1:nSampleo],:);

infoMy1 = infoM1;
infoMy2 = infoM2;
infoMy1([nSampley+1:end],:) = [];
infoMy2([nSampley+1:end],:) = [];

infoMs1 = infoM1(:,[1:nAs1]);
infoMs1([nSamples+1:end],:) = [];
infoMs2 = infoM2(:,[1:nAs2]);
infoMs2([nSamples+1:end ],:) = [];

cs1 = zeros(nAs1,nd-1);
cs2 = zeros(nAs1,nd-1);
cy1 = (infoMy1'*infoMy1)\(infoMy1'*tmpY1);
cy2 = (infoMy2'*infoMy2)\(infoMy2'*tmpY2);
for i=1:nd-1 
    cs1(:,i) = (infoMs1'*infoMs1)\(infoMs1'*tmpS1(:,i));
    cs2(:,i) = (infoMs2'*infoMs2)\(infoMs2'*tmpS2(:,i));    
end 
cs1(1,:) = 0;
cs2(1,:) = 0;

%% Sensitivity analysis 
% first moment
jcb1 = zeros(nd-1,nd-1);
jcb2 = zeros(nd-1,nd-1);
jcb1 = diag([1/dv(1),1/dv(3),1/dv(4)]);
jcb2 = diag([1/dv(2),1/dv(3),1/dv(4)]);
 
meaOfsc1 = zeros(1,nd-1);
meaOfsc2 = zeros(1,nd-1);
sen1m1 = zeros(1,nd-1);
sen1m2 = zeros(1,nd-1);

for i = 1:nd-1
    for j = 1:min(nA1,nAs1)
    sen1m1(1,i) = sen1m1(1,i) + cy1(j)*cs1(j,i);
    sen1m2(1,i) = sen1m2(1,i) + cy2(j)*cs2(j,i);
    end
end

for i =1:nd-1
    meaOfsc1(1,i) = cs1(1,i);
    meaOfsc2(1,i) = cs2(1,i);    
end 
meaOfsc1 = meaOfsc1*jcb1;
meaOfsc2 = meaOfsc2*jcb2;
sen1m1 = sen1m1*jcb1;
sen1m2 = sen1m2*jcb2;
% second moment 
sen2m1 = zeros(1,nd-1);
sen2m2 = zeros(1,nd-1);
triON1 = zeros(nA1,nA1,nAs1);
triON2 = zeros(nA2,nA2,nAs2);

for i=1:nA1 %nA=m
        for j=1:nA1 %nA=m
            for k=1:nAs1 %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m1(1,:) = sen2m1(1,:) + cy1(i)*cy1(j)*cs1(k,:);   
                    sen2m2(1,:) = sen2m2(1,:) + cy2(i)*cy2(j)*cs2(k,:);
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m1(1,:) = sen2m1(1,:) +  cy1(i)*cy1(j)*cs1(k,:);  
                                    sen2m2(1,:) = sen2m2(1,:) +  cy2(i)*cy2(j)*cs2(k,:); end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m1(1,:) = sen2m1(1,:) +  cy1(i)*cy1(j)*cs1(k,:);     
                                   sen2m2(1,:) = sen2m2(1,:) +  cy2(i)*cy2(j)*cs2(k,:);end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m1(1,:) = sen2m1(1,:) +  cy1(i)*cy1(j)*cs1(k,:); 
                                   sen2m2(1,:) = sen2m2(1,:) +  cy2(i)*cy2(j)*cs2(k,:);   end 
                end 
            else 
                % E[PsiXPsiXPsi]
                %tmpC = (infoMo'*infoMo)\(infoMo'*tmpYo);
                index = [i,j,k];
                index1 = sort(index, 'descend');
                if (triON1(index1(1), index1(2), index1(3)) == 0)
                    triON1(index1(1), index1(2), index1(3)) = sum(infoMo1(:,i).*infoMo1(:,j).*infoMo1(:,k))/nSampleo;
                else 
                end 
                
                if (triON2(index1(1), index1(2), index1(3)) == 0)
                    triON2(index1(1), index1(2), index1(3)) = sum(infoMo2(:,i).*infoMo2(:,j).*infoMo2(:,k))/nSampleo;
                else 
                end 
                                                            
                sen2m1(1,:) = sen2m1(1,:) + cy1(i)*cy1(j)*cs1(k,:)*triON1(index1(1), index1(2), index1(3));
                sen2m2(1,:) = sen2m2(1,:) + cy2(i)*cy2(j)*cs2(k,:)*triON2(index1(1), index1(2), index1(3));
            end 
        end 
    end 
end 

sen2m1 = sen2m1*jcb1;
sen2m2 = sen2m2*jcb2;

varY1 = sum(cy1(2:end).^2);
meanY1 = cy1(1);
varY2 = sum(cy2(2:end).^2);
meanY2 = cy2(1);

c(1) = meanY1 + 3*sqrt(varY1);
c(2) = meanY2 + 3*sqrt(varY2);

ceq = [];
dc1 = sen1m1' + (3/2)*(1/sqrt(varY1))*(sen2m1-2*meanY1*sen1m1)';
dc2 = sen1m2' + (3/2)*(1/sqrt(varY2))*(sen2m2-2*meanY2*sen1m2)';
index1 = [1,3,4];
index2 = [2,3,4];
DC = zeros(nd,2);
DC(index1,1) = dc1;
DC(index2,2) = dc2;

DCeq = [];
disp('inequality:')
disp(c)
save(FilNam2, 'infoM1', 'infoM2','triON1', 'triON2');
% end % End of function 
else 
% create triple product of O.G. 
%% Expansion coefficient (Y)
% y1 function output (nSampleo by 1)
load(FilNam2);
% create triple product of O.G. 
%% Expansion coefficient (Y)
% y1 function output (nSampleo by 1)
tmpY1 = zeros(nSampley,1);
tmpY2 = zeros(nSampley,1);
for L=1:nSampley
        if (L < nSampley +1) 
            tmpY1(L,1) = responY1(x(L,:), muTr);
            tmpY2(L,1) = responY2(x(L,:), muTr);
        end 
end 

tmpS1 = zeros(nSamples, nd-1);
tmpS2 = zeros(nSamples, nd-1);

for L = 1:nSample 
    if  (L <nSamples + 1) 
        scoredk1 = scoreGauss(x(L,:), mu, cov,'sub1')';
        scoredk2 = scoreGauss(x(L,:), mu, cov,'sub2')';
        for p =1:3
            tmpS1(L,p) = scoredk1(p);
            tmpS2(L,p) = scoredk2(p);
        end  
    end 
end 
% monomial basis for y1 and y2
M1 = M;
M2 = M;
M1(:,[INDEX1]) = [];
M2(:,[INDEX2]) = [];
infoM1 = (QQ2*M1')';
infoM2 = (QQ3*M2')';


infoMo1 = infoM1([1:nSampleo],:);
infoMo2  = infoM2([1:nSampleo],:);

infoMy1 = infoM1;
infoMy2 = infoM2;
infoMy1([nSampley+1:end],:) = [];
infoMy2([nSampley+1:end],:) = [];

infoMs1 = infoM1(:,[1:nAs1]);
infoMs1([nSamples+1:end],:) = [];
infoMs2 = infoM2(:,[1:nAs2]);
infoMs2([nSamples+1:end ],:) = [];

cs1 = zeros(nAs1,nd-1);
cs2 = zeros(nAs1,nd-1);
cy1 = (infoMy1'*infoMy1)\(infoMy1'*tmpY1);
cy2 = (infoMy2'*infoMy2)\(infoMy2'*tmpY2);
for i=1:nd-1 
    cs1(:,i) = (infoMs1'*infoMs1)\(infoMs1'*tmpS1(:,i));
    cs2(:,i) = (infoMs2'*infoMs2)\(infoMs2'*tmpS2(:,i));    
end 
cs1(1,:) = 0;
cs2(1,:) = 0;
%% Sensitivity analysis 
% first moment
jcb1 = zeros(nd-1,nd-1);
jcb2 = zeros(nd-1,nd-1);

jcb1 = diag([1/dv(1),1/dv(3),1/dv(4)]);
jcb2 = diag([1/dv(2),1/dv(3),1/dv(4)]);

meaOfsc1 =  zeros(1,nd-1);
meaOfsc2 =  zeros(1,nd-1);

sen1m1 = zeros(1,nd-1);
sen1m2 = zeros(1,nd-1);

for i = 1:nd-1
    for j = 1:min(nA1,nAs1)
    sen1m1(1,i) = sen1m1(1,i) + cy1(j)*cs1(j,i);
    sen1m2(1,i) = sen1m2(1,i) + cy2(j)*cs2(j,i);
    end
end
for i =1:nd-1
    meaOfsc1(1,i) = cs1(1,i);
    meaOfsc2(1,i) = cs2(1,i);    
end 
meaOfsc1 = meaOfsc1*jcb1;
meaOfsc2 = meaOfsc2*jcb2;
sen1m1 = sen1m1*jcb1;
sen1m2 = sen1m2*jcb2;
% scond moment 
sen2m1 = zeros(1,nd-1);
sen2m2 = zeros(1,nd-1);

for i=1:nA1 %nA=m
        for j=1:nA1 %nA=m
            for k=1:nAs1 %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m1(1,:) = sen2m1(1,:) + cy1(i)*cy1(j)*cs1(k,:);   
                    sen2m2(1,:) = sen2m2(1,:) + cy2(i)*cy2(j)*cs2(k,:);
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m1(1,:) = sen2m1(1,:) +  cy1(i)*cy1(j)*cs1(k,:);  
                                    sen2m2(1,:) = sen2m2(1,:) +  cy2(i)*cy2(j)*cs2(k,:); end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m1(1,:) = sen2m1(1,:) +  cy1(i)*cy1(j)*cs1(k,:);     
                                   sen2m2(1,:) = sen2m2(1,:) +  cy2(i)*cy2(j)*cs2(k,:);end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m1(1,:) = sen2m1(1,:) +  cy1(i)*cy1(j)*cs1(k,:); 
                                   sen2m2(1,:) = sen2m2(1,:) +  cy2(i)*cy2(j)*cs2(k,:);   end 
                end 
            else 
                % E[PsiXPsiXPsi]
                %tmpC = (infoMo'*infoMo)\(infoMo'*tmpYo);
                index = [i,j,k];
                index1 = sort(index, 'descend');
             
                sen2m1(1,:) = sen2m1(1,:) + cy1(i)*cy1(j)*cs1(k,:)*triON1(index1(1), index1(2), index1(3));
                sen2m2(1,:) = sen2m2(1,:) + cy2(i)*cy2(j)*cs2(k,:)*triON2(index1(1), index1(2), index1(3));
            end 
        end 
    end 
end 

sen2m1 = sen2m1*jcb1;
sen2m2 = sen2m2*jcb2;

varY1 = sum(cy1(2:end).^2);
meanY1 = cy1(1);
varY2 = sum(cy2(2:end).^2);
meanY2 = cy2(1);

c(1) = meanY1 + 3*sqrt(varY1);
c(2) = meanY2 + 3*sqrt(varY2);

ceq = [];
dc1 = sen1m1' + (3/2)*(1/sqrt(varY1))*(sen2m1-2*meanY1*sen1m1)';
dc2 = sen1m2' + (3/2)*(1/sqrt(varY2))*(sen2m2-2*meanY2*sen1m2)';
index1 = [1,3,4];
index2 = [2,3,4];
DC = zeros(nd,2);
DC(index1,1) = dc1;
DC(index2,2) = dc2;

DCeq = [];
disp('inequality:')
disp(c)

end 






            

















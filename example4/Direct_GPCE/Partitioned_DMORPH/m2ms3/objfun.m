%% ========================================================================
% Example 3: Function of objective function and its grandient (Direct GPCE w/partitioned DMORPH)   
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Input required: design variables (dv) 
%% ========================================================================
function [valObj , vObjGrad] = objfun(dv)
%clear all
%dv = ones(1,10)*30;
global cntObj stat0 statf sopt
%cntObj = 0;
% tic
double precision;
%% Initialization
% number of variables
cntObj = cntObj + 1;
N = 10;
nd = 10; % design parameter size 
m = 2; % ON degree for generic function 
ms = 3; % ON degree for score function
beta = 0.05;
% Card. of GPCE coefficients 
nA = nchoosek(N+m, m);
nAs = nchoosek(N+ms, ms);
A = max(nA, nAs);
LY = 1;
FilNam1 = sprintf('frtstm3.mat');
FilNam2 = sprintf('secstm3obj.mat');
% Define parameters of distributiosn
% normalized mean vector for N variables 
NMN = ones(1,N);
% transform to original mean vector 
MN = [dv(1), dv(2), dv(3), dv(4), dv(5), dv(6), dv(7), dv(8), dv(9), dv(10)];
% standard deviatiaon from
cov= 0.05; % coeff. of variation 
% SIG(11) is unnecessary 
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
% Covariance matrix for X1~10
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

%% First iteration 
if (cntObj == 1)
load(FilNam1); %% load ID, QQ 
load(FilNam2);
% if basis order is not enough, then increase the sampling size than five times of it. 
nSampley = 100;
nSamples = nAs*10;
nSample = max(nSampley, nSamples);
nSampleo = max(2000000, nSample);

nSample = nSampleo; 
nP = ceil(min(nA/1.5, nSampley/3));
% Sample generation for z1~10
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z = qrand(q,nSample);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);

% transformation for x1~x10
trfo = chol(cova,'lower');
x(:,1:nd) = (trfo*x(:,1:nd)')';

% lognormal 
for i=1:nd
    x(:,i) = exp(x(:,i) + mpar(i));
end 

%% Expansion coefficient (Y)
% standard least squares regression 
% information matrix

% response-values  
rsvl = zeros(nSampley,1);

% score function valued response 
rssc = zeros(nSamples,nd);

MNB = zeros(nSample, A); %monomial bases 
for i=1: A
    chkID = ID(i,:); %chkID is the same as how size of order  ex) X^2, X^3  
    nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
    nZero = length(nZeroID);
    if (nZero == 0)
        MNB(:,i) = 1;
    end 
    if (nZero == 1)
        MNB(:,i) = (x(:,nZeroID).^chkID(nZeroID));
    end 
    if (nZero == 2)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        MNB(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2));
    end 
    if (nZero == 3)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        MNB(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3));
    end 
    if (nZero == 4)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        MNB(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4));
    end 
    if (nZero == 5)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        id5 = nZeroID(5);        
        MNB(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5));
    end 
    if (nZero == 6)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        id5 = nZeroID(5);   
        id6 = nZeroID(6);           
        MNB(:,i) = (x(:,id1).^chkID(id1)).*(x(:,id2).^chkID(id2)).*(x(:,id3).^chkID(id3)).*(x(:,id4).^chkID(id4)).*(x(:,id5).^chkID(id5)).*(x(:,id6).^chkID(id6));
    end     	
end      

ONB = ORN(1:A,1:A)*MNB';    
INFM = ONB';
for L = 1:nSample 
    if (L < nSampley +1) 
        [~,~, w] = RSPSF_TRUSSf2(x(L,:), MN, 'obj');
        rsvl(L,1) = w;
    end 
    if  (L <nSamples + 1) 
        scoredk = scfcn_lgnrm(x(L,:), NMN, mpar, cova)';
        for p =1:nd
            rssc(L,p) = scoredk(p);
        end  
    end 
end

INFMy = INFM([1:nSampley],[1:nA]);
INFMs = INFM([1:nSamples],[1:nAs]);
INFMo = INFM([1:nSampleo],[1:nAs]);
initCFNy = zeros(nA, 1);

CFNs = zeros(nAs,nd);

% Obtain initial solution by SLS
INFMp = (INFMy(:,[1:nP]))'*INFMy;

initCFNy(:,1) = INFMp'*inv(INFMp*INFMp')*(INFMy(:,[1:nP])'*rsvl(:,1));
 
%% D-MORPH regression
% Generalized inverse of INFMp
[U S V] = svd(INFMp);
nRow = size(S,1);
int = 0;
for i=1:nRow 
	if (abs(S(i,i)) > 1E-7)
	int = int + 1;
	end 
end

for i=1:int
	S(i,i) = 1/S(i,i);
end 

giIMp = V*S'*U';
OPj = eye(nA)-giIMp*INFMp;
B = zeros(nA, nA);

for i=1:nA
	if (i>nP)
	B(i,i) = 1;
	end 
end

[UF SF VF] = svd(OPj*B);
cnt = 0;
E = 1;
while (E>1e-10)
	cnt = cnt + 1;
	E = abs(SF(cnt,cnt));
end
r = cnt;
wy = zeros(nA, LY);
for i=1:LY
	wy(:,i)=VF(:,r:end)*inv(UF(:,r:end)'*VF(:,r:end))*UF(:,r:end)'*initCFNy(:,i);
end

for i=nP+1:nA
    chkm1 = sum(ID(i-1,:)); % order at previous iter.  
    chkm2 = sum(ID(i,:)); % order at current iter. 
    
    if (i == (nP + 1))
        if (chkm1 == chkm2)
            for j=1:LY
                crt = nchoosek(N+chkm1-1, chkm1-1);			
                wy(i,j) = abs(mean((wy(crt+1:nP))));
            end 
        else 
            for j=1:LY
                crt = nchoosek(N+chkm1-1, chkm1-1);			
                wy(i,j) = beta*abs(mean((wy(crt+1:nP))));
            end
        end 
    else 
        if (chkm1 == chkm2)
            for j=1:LY
                wy(i,j) = wy(i-1,j);
            end 
        else 
            for j=1:LY
                wy(i,j) = beta*(wy(i-1,j));
            end 
        end 
    end 
end 

[UF SF VF] = svd(OPj);
% Find r
cnt = 0;
E = 1;
while (E > 1e-6)
    cnt = cnt + 1;
    E = abs(SF(cnt,cnt));
end 
r = cnt;
SFr = zeros(r-1,r-1);
for i=1:r-1
    SFr(i,i) = 1/SF(i,i);
end 

for i=1:1
    CFNy(:,i) = VF(:,r:end)*inv(UF(:,r:end)'*VF(:,r:end))*UF(:,r:end)'*initCFNy(:,i)...
    +VF(:,1:r-1)*inv(UF(:,1:r-1)'*VF(:,1:r-1))*SFr*UF(:,1:r-1)'*OPj*wy(:,i);
end 

for i=1:nd 
    CFNs(:,i) = (INFMs'*INFMs)\(INFMs'*rssc(:,i));
end

CFNs(1,:) = 0;
%% Sensitivity analysis
% first moment
jcb = zeros(nd,nd);
for i = 1:nd
    jcb(i,i) = 1/dv(i);
end 

sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) + CFNy(j)*CFNs(j,i);
    end
end 

sen1m = sen1m*jcb;
% second moment 
sen2m = zeros(1,nd);
TRON = zeros(nA, nA, nAs);

for i=1:nA %nA=m
    disp(i)
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(1,:) = sen2m(1,:) + CFNy(i)*CFNy(j)*CFNs(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  CFNy(i)*CFNy(j)*CFNs(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  CFNy(i)*CFNy(j)*CFNs(k,:);  end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  CFNy(i)*CFNy(j)*CFNs(k,:);  end 
                end 
                else 
                % E[PsiXPsiXPsi]
                index = [i,j,k];
				index1 = sort(index, 'ascend');
				if (TRON(index1(1), index1(2), index1(3)) == 0)
				TRON(index1(1), index1(2), index1(3)) = sum(INFMo(:,i).*INFMo(:,j).*INFMo(:,k))/nSampleo;
                else 
				end 
                %tmpC = (infoMo'*infoMo)\(infoMo'*tmpYo);
                sen2m(1,:) = sen2m(1,:) + CFNy(i)*CFNy(j)*CFNs(k,:)*TRON(index1(1), index1(2), index1(3));
            end 
        end 
    end 
end 

sen2m = sen2m*jcb;

varaOfobj = sum(CFNy(2:end).^2);
meaOfobj = CFNy(1);

w1 =12589.4022163097; 
w2 = 334.41;
valObj = 0.5*meaOfobj/w1 + 0.5*sqrt(varaOfobj)/w2;

vObjGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varaOfobj))*(sen2m - 2*meaOfobj*sen1m ))';
disp(dv) 
disp('var of obj:');
disp(varaOfobj)
disp('mean of obj:');
disp(meaOfobj)

% data saving for 2nd step 
save(FilNam2, 'INFM', 'x', 'TRON', 'CFNs', 'w1', '-v7.3');

% from 2nd iterations 
else %(cntObj~=1)
    load(FilNam1);
    load(FilNam2);
    nSampley = 100;
    nSamples = nAs*10;
    nSample = max(nSampley, nSamples);
    nSampleo = max(2000000, nSample);
    nSample = nSampleo;
nP = ceil(min(nA/1.5, nSampley/3));
    % (Part needed for updating at next steps) 
rsvl = zeros(nSampley,1);
for L = 1:nSample 
    if (L < nSampley +1) 
        [~,~, w] = RSPSF_TRUSSf2(x(L,:), MN, 'obj');
        rsvl(L,1) = w;
    end 
end 
INFMy = INFM([1:nSampley],[1:nA]);
INFMs = INFM([1:nSamples],[1:nAs]);
INFMo = INFM([1:nSampleo],[1:nA]);

CFNy = zeros(nA, 1);
INFMp = (INFMy(:,[1:nP]))'*INFMy; %nP X nA
for i=1:1
	initCFNy(:,i) = INFMp'*inv(INFMp*INFMp')*(INFMy(:,[1:nP])'*rsvl(:,i));
end 
%% D-MORPH regression 
% Generalized inverse of INFMp
[U S V] = svd(INFMp);
nRow = size(S,1);
int = 0;
for i=1:nRow 
	if (abs(S(i,i)) > 1E-7)
	int = int + 1;
	end 
end

for i=1:int
	S(i,i) = 1/S(i,i);
end 

giIMp = V*S'*U';
OPj = eye(nA)-giIMp*INFMp;
B = zeros(nA, nA);

for i=1:nA
	if (i>nP)
	B(i,i) = 1;
	end 
end

[UF SF VF] = svd(OPj*B);
cnt = 0;
E = 1;
while (E>1e-10)
	cnt = cnt + 1;
	E = abs(SF(cnt,cnt));
end
r = cnt;
wy = zeros(nA, LY);
for i=1:LY
	wy(:,i)=VF(:,r:end)*inv(UF(:,r:end)'*VF(:,r:end))*UF(:,r:end)'*initCFNy(:,i);
end

for i=nP+1:nA
    chkm1 = sum(ID(i-1,:)); % order at previous iter.  
    chkm2 = sum(ID(i,:)); % order at current iter. 
    
    if (i == (nP + 1))
        if (chkm1 == chkm2)
            for j=1:LY
                crt = nchoosek(N+chkm1-1, chkm1-1);			
                wy(i,j) = abs(mean((wy(crt+1:nP))));
            end 
        else 
            for j=1:LY
                crt = nchoosek(N+chkm1-1, chkm1-1);			
                wy(i,j) = beta*abs(mean((wy(crt+1:nP))));
            end
        end 
    else 
        if (chkm1 == chkm2)
            for j=1:LY
                wy(i,j) = wy(i-1,j);
            end 
        else 
            for j=1:LY
                wy(i,j) = beta*(wy(i-1,j));
            end 
        end 
    end 
end 

[UF SF VF] = svd(OPj);
% Find r
cnt = 0;
E = 1;
while (E > 1e-6)
    cnt = cnt + 1;
    E = abs(SF(cnt,cnt));
end 
r = cnt;
SFr = zeros(r-1,r-1);
for i=1:r-1
    SFr(i,i) = 1/SF(i,i);
end 

for i=1:1
    CFNy(:,i) = VF(:,r:end)*inv(UF(:,r:end)'*VF(:,r:end))*UF(:,r:end)'*initCFNy(:,i)...
    +VF(:,1:r-1)*inv(UF(:,1:r-1)'*VF(:,1:r-1))*SFr*UF(:,1:r-1)'*OPj*wy(:,i);
end 

 %% Sensitivity analysis 
% first moment
jcb = zeros(nd,nd);
for i = 1:nd
    jcb(i,i) = 1/dv(i);
end 

sen1m = zeros(1,nd);
for i = 1:nd
    for j = 1:min(nA,nAs)
    sen1m(1,i) = sen1m(1,i) + CFNy(j)*CFNs(j,i);
    end
end 

sen1m = sen1m*jcb;

% second moment 
sen2m = zeros(1,nd);
    for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(1,:) = sen2m(1,:) + CFNy(i)*CFNy(j)*CFNs(k,:);                        
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(1,:) = sen2m(1,:) +  CFNy(i)*CFNy(j)*CFNs(k,:);  end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(1,:) = sen2m(1,:) +  CFNy(i)*CFNy(j)*CFNs(k,:);  end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(1,:) = sen2m(1,:) +  CFNy(i)*CFNy(j)*CFNs(k,:);  end 
                end 
                else  
                % E[PsiXPsiXPsi]
                index = [i,j,k];
				index1 = sort(index, 'ascend');
                sen2m(1,:) = sen2m(1,:) + CFNy(i)*CFNy(j)*CFNs(k,:)*TRON(index1(1), index1(2), index1(3));
            end 
        end 
    end 
    end     

sen2m = sen2m*jcb;
    
varaOfobj = sum(CFNy(2:end).^2);
meaOfobj = CFNy(1);

w1 =12589.4022163097; 
w2 = 334.41;
valObj = 0.5*meaOfobj/w1 + 0.5*sqrt(varaOfobj)/w2;
%valObj = meaOfobj/w1;
%vObjGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varaOfobj))*(sen2m - 2*meaOfobj*sen1m))'; 
vObjGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varaOfobj))*(sen2m - 2*meaOfobj*sen1m ))'; 
disp(dv) 
disp('var of obj :');
disp(varaOfobj);
disp('mean of obj:');
disp(meaOfobj);

end % end of function 

% record stat. info. at the initial and final design. 
switch sopt
    case 'pre'
        stat0 = [meaOfobj, sqrt(varaOfobj)];
    case 'post'
        statf = [meaOfobj, sqrt(varaOfobj)];
    otherwise 
end 



            

















%% ========================================================================
% Example 4: Function of inequality constraint function and its grandient (Direct GPCE w/DMORPH)   
% Input: design variables (dv) 
% Output: constraint values and design sensitivities
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ======================================================================== 
function[cnst, ceq, grdnV,  grdneq] = confun(dv)

global cntCon cnst1 diffcnst
double precision;
%% Initialization
cntCon = cntCon + 1;
N = 10; % number of random variables 
nd = 10; % design parameter size 
m = 2; % order of GPCE for generic function 
ms = 3; % order of GPCE for score function
beta = 0.05;
% Card. of GPCE coefficients 
nA = nchoosek(N+m, m);
nAs = nchoosek(N+ms, ms);
A = max(nA,nAs);

FilNam1 = sprintf('frtstm3.mat');
FilNam2 = sprintf('secstm3con.mat');

% normalized mean vector for N variables 
NMN = ones(1,N);
% transform to original mean vector 
MN = [dv(1), dv(2), dv(3), dv(4), dv(5), dv(6), dv(7), dv(8), dv(9), dv(10)];

% standard deviatiaon from
cov= 0.05; % coeff. of variation 
SIG = [cov, cov, cov, cov, cov, cov, cov, cov, cov, cov];

% log-normal parameters
for i = 1:nd
mpar(i) = log(NMN(i)^2/sqrt(SIG(i)^2 + NMN(i)^2)); % paramter mu 
spar(i) = sqrt(log(1 + ((SIG(i))^2)/(NMN(i)^2))); % paramter sig
end

% Correlation coefficient matrix in normal dist. 
cor = zeros(nd,nd);
for i=1:nd
    for j=i:nd
        if (i~=j)
            cor(i,j) = 0.2;
        end 
        if (i==j)
            cor(i,j) = 1;
        end 
    end 
end 
for i=1:nd
    for j=i:nd
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

%% First iteration
if (cntCon == 1)
load(FilNam1);   
load(FilNam2); 
clear rssc;
clear x;
nSampley = 100;
nSamples = nAs*10;
nSample = max(nSampley, nSamples); 
nSampleo = max(2000000, nSample);

nSample = nSampleo; 
nP = ceil(min(nA/1.5, nSampley/3));
% Sampling for Z1~10
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z = qrand(q,nSample);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);

% transformation z1~z10 to x1~x10
trfo = chol(cova,'lower');
x(:,1:nd) = (trfo*x(:,1:nd)')';

% lognormal 
for i=1:nd
    x(:,i) = exp(x(:,i) + mpar(i));
end 
% Monomial basis (MNB) 
MNB = zeros(nSample, A);
for i=1: A
    chkID = ID(i,:); %chkID: ex) (x1^(2), x2^(3))->(2,3)  
    nZeroID = find(chkID~=0); 
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

cnt = 1;

for L = 1:nSample
    if (L < nSampley + 1) 
        [Y, LY, ~] = RSPSF_TRUSSf2(x(L,:), MN, 'cnstn');
		if (cnt == 1), rsvl = zeros(nSampley, LY); end
        cnt = cnt + 1;
        % 'cnstn' option create two performance function values of truss (check inside) 
        rsvl(L,:) = Y;
       %rsvl2(L,1) = Y2;
    end 
    if  (L <nSamples + 1) 
        scoredk = scfcn_lgnrm(x(L,:), NMN, mpar, cova)';
        for p =1:nd
            rssc(L,p) = scoredk(p);
        end  
    end 
end 

INFMy = INFM(1:nSampley, 1:nA);
INFMs = INFM(1:nSamples, 1:nAs);
INFMo = INFM(1:nSampleo, 1:A);
initCFNy = zeros(nA, LY);

CFNs = zeros(nAs,nd);
CFNy = zeros(nA, LY);
% Obtain initial solution by SLS
INFMp = (INFMy(:,[1:nP]))'*INFMy; %nP X nA
for i=1:LY
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
    chkm1 = sum(ID(i-1,:)); 
    chkm2 = sum(ID(i,:));
    
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

for i=1:LY
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

sen1m = zeros(LY,nd);
for k=1:LY
	for i = 1:nd
		for j = 1:min(nA,nAs)
		sen1m(k,i) = sen1m(k,i) + CFNy(j,k)*CFNs(j,i);
		end
	end
end 
 
sen1m = sen1m*jcb;
% second moment 
sen2m = zeros(LY,nd);
TRON = zeros(nA, nA, nAs);
for p = 1:LY
    for i=1:nA %nA=m
        if (p==1), disp(i); end 
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:);   
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(p,:) = sen2m(p,:) +  CFNy(i,p)*CFNy(j,p)*CFNs(k,:); end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(p,:) = sen2m(p,:) +  CFNy(i,p)*CFNy(j,p)*CFNs(k,:); end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:);  end 
                end 
            else 
                index = [i,j,k];
                if (nA <= nAs), index1 = sort(index, 'ascend'); end 
                if (nA > nAs), index1 = sort(index, 'descend'); end
                if (TRON(index1(1), index1(2), index1(3)) == 0)
				TRON(index1(1), index1(2), index1(3)) = sum(INFMo(:,i).*INFMo(:,j).*INFMo(:,k))/nSampleo; end 
                sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:)*TRON(index1(1), index1(2), index1(3));
            end 
        end 
    end 
    end 
end 

sen2m = sen2m*jcb;

varaOfcnst  = sum(CFNy(2:end,:).^2);
meaOfcnst = CFNy(1,:);

cnst = (meaOfcnst + 3.*sqrt(varaOfcnst));
cnst1 = cnst(:);
ceq = [];
for i=1:LY
	grdnV(:,i) = sen1m(i,:)' + (3/2)*(1/sqrt(varaOfcnst(i)))*(sen2m(i,:)-2*meaOfcnst(i)*sen1m(i,:) )';
end 
diffcnst = grdnV(:);
grdneq = [];
disp('Inequality:');
disp(cnst);
% end % End of function 
save(FilNam2, 'TRON', 'INFM', 'x', 'CFNs', '-v7.3');


else %% 2nd or higher iterations 

load(FilNam1);   
load(FilNam2); 

nSampley = 100;
nSamples = nAs*10;
nSample = max(nSampley, nSamples); 
nSampleo = max(2000000, nSample);

nSample = nSampleo; 
nP = ceil(min(nA/1.5, nSampley/3));

cnt = 1;

for L = 1:nSample
    if (L < nSampley + 1) 
        [Y, LY, ~] = RSPSF_TRUSSf2(x(L,:), MN, 'cnstn');
		if (cnt == 1), rsvl = zeros(nSampley, LY); end
        cnt = cnt + 1;
        % 'cnstn' option create two performance function values of truss (check inside) 
        rsvl(L,:) = Y;
    end 
end 

INFMy = INFM(1:nSampley, 1:nA);
INFMs = INFM(1:nSamples, 1:nAs);
INFMo = INFM(1:nSampleo, 1:A);
initCFNy = zeros(nA, LY);

CFNy = zeros(nA, LY);
% Obtain initial solution by SLS
INFMp = (INFMy(:,[1:nP]))'*INFMy; %nP X nA
for i=1:LY
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
    chkm1 = sum(ID(i-1,:)); 
    chkm2 = sum(ID(i,:));
    
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

for i=1:LY
    CFNy(:,i) = VF(:,r:end)*inv(UF(:,r:end)'*VF(:,r:end))*UF(:,r:end)'*initCFNy(:,i)...
    +VF(:,1:r-1)*inv(UF(:,1:r-1)'*VF(:,1:r-1))*SFr*UF(:,1:r-1)'*OPj*wy(:,i);
end 
        
%% Sensitivity analysis 
% first moment
jcb = zeros(nd,nd);
for i = 1:nd
    jcb(i,i) = 1/dv(i);
end 

sen1m = zeros(LY,nd);
for k=1:LY
	for i = 1:nd
		for j = 1:min(nA,nAs)
		sen1m(k,i) = sen1m(k,i) + CFNy(j,k)*CFNs(j,i);
		end
	end
end 
 
sen1m = sen1m*jcb;
% second moment 
sen2m = zeros(LY,nd);
for p = 1:LY
    for i=1:nA %nA=m
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1))
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:);   
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(p,:) = sen2m(p,:) +  CFNy(i,p)*CFNy(j,p)*CFNs(k,:); end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(p,:) = sen2m(p,:) +  CFNy(i,p)*CFNy(j,p)*CFNs(k,:); end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:);  end 
                end 
            else 
                index = [i,j,k];
                if (nA <= nAs), index1 = sort(index, 'ascend'); end 
                if (nA > nAs), index1 = sort(index, 'descend'); end
                sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:)*TRON(index1(1), index1(2), index1(3));
            end 
        end 
    end 
    end 
end 

sen2m = sen2m*jcb;

varaOfcnst  = sum(CFNy(2:end,:).^2);
meaOfcnst = CFNy(1,:);

cnst = (meaOfcnst + 3.*sqrt(varaOfcnst));
cnst1 = cnst(:);
ceq = [];
for i=1:LY
	grdnV(:,i) = sen1m(i,:)' + (3/2)*(1/sqrt(varaOfcnst(i)))*(sen2m(i,:)-2*meaOfcnst(i)*sen1m(i,:) )';
end 
diffcnst = grdnV(:);
grdneq = [];
disp('Inequality:');
disp(cnst);
% End of function 
end 





            

















%% ========================================================================
% Example 4: Function of inequality constraint function and its grandient (Single-step GPCE)   
% Input: design variables (dv) 
% Output: constraint value and design sensitivities 
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ======================================================================== 
function[cnst, ceq, grdnV,  grdneq] = sconfun(dv)

global cntCon cnst1 diffcnst jj
double precision;
%% Initialization
% number of variables
cntCon = cntCon + 1;
jj= jj + 1;
N = 10; % number of random variables 
nd = 10; % number of design variables 
m = 1; % order of GPCE for performance function 
ms = 3; % order of GPCE for score function

nA = nchoosek(N+m, m); % number of GPCE coefficients for performance function 
nAs = nchoosek(N+ms, ms); % number of GPCE coefficients for score function 
A = max(nA,nAs); % maximum number of GPCE coefficients 

FilNam1 = sprintf('frtstm3cor2.mat');
FilNam2 = sprintf('secstm3con.mat');


NMN = ones(1,N); % normalized mean vector 
% 
MN = [dv(1), dv(2), dv(3), dv(4), dv(5), dv(6), dv(7), dv(8), dv(9), dv(10)]; % original mean vector

cov= 0.05; % coeff. of variation 
SIG = [cov, cov, cov, cov, cov, cov, cov, cov, cov, cov];

%% Set log-normal parameters
for i = 1:nd
mpar(i) = log(NMN(i)^2/sqrt(SIG(i)^2 + NMN(i)^2)); % paramter mu 
spar(i) = sqrt(log(1 + ((SIG(i))^2)/(NMN(i)^2))); % paramter sig
end

%% Set correlation coefficient matrix in normal dist. 
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

% Construct covariance matrix
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
if (jj ~= 1)
load(FilNam2); end    
nSampley = nA*3; % sampling size for performance func. 
nSamples = nAs*10; % sampling size for score func. 
nSample = max(nSampley, nSamples); 
nSampleo = max(2000000, nSample);

nSample = nSampleo; 

%% Construct QMCS 
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z = qrand(q,nSample);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);

%% transform to correlated random vector 
trfo = chol(cova,'lower');
x(:,1:nd) = (trfo*x(:,1:nd)')';

% transform to lognormal random vector 
for i=1:nd
    x(:,i) = exp(x(:,i) + mpar(i));
end 

%% Expansion coefficient (Y)
% y function output (nSampleo by 1)
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

ONB = ORN(1:A,1:A)*MNB'; % construct ON basis 
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
	if (jj == 1)
    if  (L <nSamples + 1) 
        scoredk = scfcn_lgnrm(x(L,:), NMN, mpar, cova)';
        for p =1:nd
            rssc(L,p) = scoredk(p);
        end  
    end, end  
end 

INFMy = INFM(1:nSampley, 1:nA);
INFMs = INFM(1:nSamples, 1:nAs);
INFMo = INFM(1:nSampleo, 1:A);

if (jj== 1)
CFNs = zeros(nAs,nd);
for i=1:nd 
    CFNs(:,i) = (INFMs'*INFMs)\(INFMs'*rssc(:,i));
end 
CFNs(1,:) = 0;
end 

CFNy = zeros(nA, LY);

for i=1:LY
	CFNy(:,i) = (INFMy'*INFMy)\(INFMy'*rsvl(:,i));	
end 



%% Sensitivity analysis 

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
if (jj == 1)
TRON = zeros(nA, nA, nAs); end 
for p = 1:LY
    for i=1:nA %nA=m
        if (p==1), disp(i); end 
        for j=1:nA %nA=m
            for k=1:nAs %nA=m'
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
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
                % E[PsiXPsiXPsi]
                %tmpC = (infoMo'*infoMo)\(infoMo'*tmpYo);
                index = [i,j,k];
                if (nA <= nAs), index1 = sort(index, 'ascend'); end 
                if (nA > nAs), index1 = sort(index, 'descend'); end
                if (TRON(index1(1), index1(2), index1(3)) == 0)
				TRON(index1(1), index1(2), index1(3)) = sum(INFMo(:,i).*INFMo(:,j).*INFMo(:,k))/nSampleo; end 
                %tmpC = (infoMo'*infoMo)\(infoMo'*tmpYo);
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
CFNy0 = CFNy;
dv0 = dv;

disp('Inequality:');
disp(cnst);
% end % End of function 
save(FilNam2, 'TRON', 'INFM', 'x', 'CFNs', 'CFNy0', 'dv0', 'LY', '-v7.3');

% 2nd or higher iterations 
else 
load(FilNam1);   
load(FilNam2); 

nSampley = nA*3;
nSamples = nAs*10;
nSample = max(nSampley, nSamples); 
nSampleo = max(2000000, nSample);

nSample = nSampleo; 

MNB = zeros(nSample, A); %monomial bases 
for i=1: A
    chkID = ID(i,:); %chkID is the same as how size of order  ex) X^2, X^3  
    nZeroID = find(chkID~=0); %nZeroID is the same as which of variables ex) X1, X2 
    nZero = length(nZeroID);
    if (nZero == 0)
        MNB(:,i) = 1;
    end 
    if (nZero == 1)
        MNB(:,i) = (((dv(nZeroID)/dv0(nZeroID))*x(:,nZeroID)).^chkID(nZeroID));
    end 
    if (nZero == 2)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        MNB(:,i) = (((dv(id1)/dv0(id1))*x(:,id1)).^chkID(id1)).*(((dv(id2)/dv0(id2))*x(:,id2)).^chkID(id2));
    end 
    if (nZero == 3)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        MNB(:,i) = (((dv(id1)/dv0(id1))*x(:,id1)).^chkID(id1)).*(((dv(id2)/dv0(id2))*x(:,id2)).^chkID(id2)).*(((dv(id3)/dv0(id3))*x(:,id3)).^chkID(id3));
    end 
    if (nZero == 4)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        MNB(:,i) = (((dv(id1)/dv0(id1))*x(:,id1)).^chkID(id1)).*(((dv(id2)/dv0(id2))*x(:,id2)).^chkID(id2)).*(((dv(id3)/dv0(id3))*x(:,id3)).^chkID(id3)).*(((dv(id4)/dv0(id4))*x(:,id4)).^chkID(id4));
    end 
    if (nZero == 5)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        id5 = nZeroID(5);        
        MNB(:,i) = (((dv(id1)/dv0(id1))*x(:,id1)).^chkID(id1)).*(((dv(id2)/dv0(id2))*x(:,id2)).^chkID(id2)).*(((dv(id3)/dv0(id3))*x(:,id3)).^chkID(id3)).*(((dv(id4)/dv0(id4))*x(:,id4)).^chkID(id4)).*(((dv(id5)/dv0(id5))*x(:,id5)).^chkID(id5));
    end 
    if (nZero == 6)
        id1 = nZeroID(1);
        id2 = nZeroID(2);
        id3 = nZeroID(3);
        id4 = nZeroID(4);
        id5 = nZeroID(5);   
        id6 = nZeroID(6);           
        MNB(:,i) = (((dv(id1)/dv0(id1))*x(:,id1)).^chkID(id1)).*(((dv(id2)/dv0(id2))*x(:,id2)).^chkID(id2)).*(((dv(id3)/dv0(id3))*x(:,id3)).^chkID(id3)).*(((dv(id4)/dv0(id4))*x(:,id4)).^chkID(id4)).*(((dv(id5)/dv0(id5))*x(:,id5)).^chkID(id5)).*(((dv(id6)/dv0(id6))*x(:,id6)).^chkID(id6));
    end     	
end    

cnt = 1;
LY = 11;
for k=1:LY
	for L = 1:nSample
		if (L < nSampley + 1) 
			Y = CFNy0(:,k)'*(ORN(1:nA,1:nA)*MNB(L,1:nA)');
        % 'cnstn' option create two performance function values of truss (check inside) 
			rsvl(L,k) = Y;
       %rsvl2(L,1) = Y2;
		end
	end 		
end 

INFMy = INFM(1:nSampley, 1:nA);
INFMs = INFM(1:nSamples, 1:nAs);
INFMo = INFM(1:nSampleo, 1:A);

CFNy = zeros(nA, LY);

for i=1:LY
	CFNy(:,i) = (INFMy'*INFMy)\(INFMy'*rsvl(:,i));	
end 

%% Sensitivity analysis 

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
                if ( (i == 1) || (j == 1) || (k ==1)) % table 
                if ((i==j) && (i==k) && (j==k)) 
                    sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:);   
                elseif ((i==1) && (j~=1))
                    if (j==k),  sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:); end 
                elseif ((j==1) && (i~=1))
                    if (i==k),  sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:); end 
                elseif ((k==1) && (i~=1))
                    if (i==j),  sen2m(p,:) = sen2m(p,:) + CFNy(i,p)*CFNy(j,p)*CFNs(k,:);  end 
                end 
            else 
                % E[PsiXPsiXPsi]
                %tmpC = (infoMo'*infoMo)\(infoMo'*tmpYo);
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
% end % End of function 
end 





            

















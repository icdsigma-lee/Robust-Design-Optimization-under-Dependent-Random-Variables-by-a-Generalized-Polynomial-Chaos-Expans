%% ========================================================================
% Example 3: Function of inequality constraint function and its grandient (Direct GPCE w/SLS)   
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Input required: design variables (dv) 
%% ======================================================================== 
function[cnst, ceq, grdnV,  grdneq] = confun(dv)

global cntCon
double precision;
%% Initialization
% number of variables
%cntCon = 0;
cntCon = cntCon + 1;
N = 10;
nd = 10; % design parameter size 
m = 2; % ON degree for generic function 
ms = 3; % ON degree for score function

% Card. of GPCE coefficients 
nA = nchoosek(N+m, m);
nAs = nchoosek(N+ms, ms);

FilNam1 = sprintf('frtstm3cor2.mat');
FilNam2 = sprintf('secstm3.mat');

% normalized mean vector for N variables 
NMN = ones(1,N);
% transform to original mean vector 
MN = [dv(1), dv(2), dv(3), dv(4), dv(5), dv(6), dv(7), dv(8), dv(9), dv(10)];

% standard deviatiaon from
cov= 0.05; % coeff. of variation 
SIG = [cov, cov, cov, cov, cov, cov, cov, cov, cov, cov];

% log-normal  
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

% load 
load(FilNam1);   
load(FilNam2);

nSampley = nA*3;
nSamples = nAs*10;
nSampleo = 2000000;

%nSampleo = nAo*nt;
glSample = [nSampley, nSamples, nSampleo];
nSample = max(glSample); 

%% Expansion coefficient (Y)
% y function output (nSampleo by 1)
cnt = 1;
for L=1:nSampley 
        [Y, LY, ~] = RSPSF_TRUSSf2(x(L,:), MN, 'cnstn');
		if (cnt == 1), rsvl = zeros(nSampley, LY); end
        cnt = cnt + 1;
        % 'cnstn' option create two performance function values of truss (check inside) 
        rsvl(L,:) = Y;
       %rsvl2(L,1) = Y2;
end 

INFMy = INFM([1:nSampley], [1:nA]);
INFMs = INFM([1:nSamples], [1:nAs]);
INFMo = INFM([1:nSampleo], [1:nAs]);

CFNy = zeros(nA, LY);

for i=1:LY
	CFNy(:,i) = (INFMy'*INFMy)\(INFMy'*rsvl(:,i));	
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
				index1 = sort(index, 'ascend');
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

ceq = [];
for i=1:LY
	grdnV(:,i) = sen1m(i,:)' + (3/2)*(1/sqrt(varaOfcnst(i)))*(sen2m(i,:)-2*meaOfcnst(i)*sen1m(i,:) )';
end 

grdneq = [];
disp('Inequality:');
disp(cnst);
% end % End of function 










            

















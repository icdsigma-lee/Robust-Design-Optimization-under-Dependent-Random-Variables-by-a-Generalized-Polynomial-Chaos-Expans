%% ========================================================================
% Example 3: Function of objective function and its grandient (Signle-step GPCE)   
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
% Input required: design variables (dv) 
%% ========================================================================
function [valObj , vObjGrad] = sobjfun(dv)
%clear all
global cntObj stat0 statf sopt ii
%cntObj = 0;
% tic
double precision;
%% Initialization
% number of variables
cntObj = cntObj + 1;
ii = ii + 1;
N = 10; % number of random variables 
nd = 10; % number of design variables  
m = 1; % ON degree for generic function 
ms = 3; % ON degree for score function
 
nA = nchoosek(N+m, m); % number of GPCE coefficients in performance func.
nAs = nchoosek(N+ms, ms); % number of GPCE coefficients in score func.
A = max(nA, nAs); % maximum number of GPCE coefficients 

FilNam1 = sprintf('frtstm3cor2.mat');
FilNam2 = sprintf('secstm3obj.mat');
%% Define parameters of distribution
  
NMN = ones(1,N); % normalized mean vector 

MN = [dv(1), dv(2), dv(3), dv(4), dv(5), dv(6), dv(7), dv(8), dv(9), dv(10)]; % original mean vector 

cov= 0.05; % coefficient of variation 

SIG = [cov, cov, cov, cov, cov, cov, cov, cov, cov, cov]; 

%% log-normal parameters 
for i = 1:nd
mpar(i) = log(NMN(i)^2/sqrt(SIG(i)^2 + NMN(i)^2)); % paramter mu 
spar(i) = sqrt(log(1 + ((SIG(i))^2)/(NMN(i)^2))); % paramter sig
end
 
%% Set correlation coefficient matrix in normal dist. 
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
%% Set covariance matrix for X1~10
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
load(FilNam1); %% load ID, ORN  
if (ii ~= 1), load(FilNam2); end  
% if basis order is not enough, then increase the sampling size than five times of it. 
nSampley = nA*3; % sampling size for a performance function
nSamples = nAs*10; % sampling size efor score function 
nSample = max(nSampley, nSamples);  
nSampleo = max(2000000, nSample);
nSample = nSampleo; 

%% Construct input-output data for Z1~10
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z = qrand(q,nSample);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);

%% Transform to correlated random vector 
trfo = chol(cova,'lower');
x(:,1:nd) = (trfo*x(:,1:nd)')';

% Trannsform to lognormal random vector  
for i=1:nd
    x(:,i) = exp(x(:,i) + mpar(i));
end 

%% Expansion coefficient (Y)
% standard least squares regression 
% information matrix

% output set  
rsvl = zeros(nSampley,1);

% score function valued output set 
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

%% Construct ON 
ONB = ORN(1:A,1:A)*MNB';    
INFM = ONB';

%% Construct output samples 
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

INFMy = INFM(1:nSampley,1:nA);
INFMs = INFM(1:nSamples,1:nAs);
INFMo = INFM(1:nSampleo,1:A);


CFNy = (INFMy'*INFMy)\(INFMy'*rsvl); % GPCE coefficients for a performance function
if (ii == 1)
CFNs = zeros(nAs,nd);
    for i=1:nd 
        CFNs(:,i) = (INFMs'*INFMs)\(INFMs'*rssc(:,i));  % GPCE coefficients for a score function 
    end
CFNs(1,:) = 0;
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
if (ii== 1)
TRON = zeros(nA, nA, nAs);
else 
end 

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
                if (nA <= nAs), index1 = sort(index, 'ascend'); end 
                if (nA > nAs), index1 = sort(index, 'descend'); end 
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
vObjGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varaOfobj))*(sen2m - 2*meaOfobj*sen1m))'; 

disp(dv) 
disp('var of obj:');
disp(varaOfobj)
disp('mean of obj:');
disp(meaOfobj)
CFNy0 = CFNy;
dv0 = dv;
% data saving for 2nd step 
save(FilNam2, 'INFM', 'x', 'TRON', 'CFNs', 'CFNy0', 'dv0', '-v7.3');

% 2nd or higher iterations 
else %(cntObj~=1)
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

% (Part needed for updating at next steps) 
rsvl = zeros(nSampley,1);
for L = 1:nSample 
    if (L < nSampley +1) 
        %[~,~, w] = RSPSF_TRUSSf2(x(L,:), MN, 'obj');
		w = CFNy0'*(ORN(1:nA,1:nA)*MNB(L,1:nA)');
        rsvl(L,1) = w;
    end 
end 
INFMy = INFM([1:nSampley],[1:nA]);
INFMs = INFM([1:nSamples],[1:nAs]);
INFMo = INFM([1:nSampleo],[1:A]);

CFNy = (INFMy'*INFMy)\(INFMy'*rsvl);

 %% Sensitivity analysis 

jcb = zeros(nd,nd);
for i = 1:nd
    jcb(i,i) = 1/dv(i);
end 

% first moment
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
                if (nA <= nAs), index1 = sort(index, 'ascend'); end 
                if (nA > nAs), index1 = sort(index, 'descend'); end 
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
vObjGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varaOfobj))*(sen2m - 2*meaOfobj*sen1m))'; 

disp(dv) 
disp('var of obj :');
disp(varaOfobj);
disp('mean of obj:');
disp(meaOfobj);

end % end of function 

% record stat. info. at the initial and final design. 
switch sopt
    case 'pre'
        stat0 = [meaOfobj, varaOfobj];
    case 'post'
        statf = [meaOfobj, varaOfobj];
    otherwise 
end 



            

















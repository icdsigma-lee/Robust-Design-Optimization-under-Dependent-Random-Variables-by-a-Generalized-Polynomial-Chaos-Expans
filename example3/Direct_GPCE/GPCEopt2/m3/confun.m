%% ========================================================================
% Example 3: Function of inequality constraint function and its grandient (Direct GPCE option2)   
% Input: design variables (dv) 
% Output: constraint values and design sensitivities
% written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function[cnst, ceq, grdnV,  grdneq] = confun(dv)
%dv = ones(1,10)*30;
global cntCon diffcnst cnst1
double precision;
%% Initialization
% number of variables
%cntCon = 0;
cntCon = cntCon + 1;
N = 7; % number of random variables 
nd = 4; % design parameter size 
m = 3; % order of GPCE for generic function 
% Card. of GPCE coefficients 
nA = nchoosek(N+m, m);
nA1 = nchoosek(N-2+m, m);
FilNam1 = sprintf('frtst.mat');
FilNam2 = sprintf('secstcon.mat');

% normalized mean vector for N variables 
tmu5 = 10000; %kgm^3
tmu6 = 2050000000; %Pa
tmu7 = 200000; %N
NMN = ones(1,N);
% transform to original mean vector 
MN = [dv(1), dv(2), dv(3), dv(4), tmu5, tmu6, tmu7];

% standard deviatiaon from
SIG = [0.02, 0.02, 0.02, 0.02, 0.3, 25/105, 0.25];

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
 
%% First iteration 
if (cntCon == 1)
load(FilNam1); %% load ID, QQ 
% if basis order is not enough, then increase the sampling size more than five times of it. 
nSample = nA*15;

% QMC Sampling for Z1~10
rng(123457);
p = sobolset(N,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
q = qrandstream(p);
z = qrand(q,nSample);
z1 = unifcdf(z,0,1);
x = norminv(z1,0,1);

% Covariance matrix for x1~x2 
cov12 = zeros(2,2);
for i = 1:2
    for j = i:2
        cov12(i,j) = cor(i,j)*SIG(i)*SIG(j);
    end 
end 
for i = 1:2
    for j=i:2
        cov12(j,i) = cov12(i,j);
    end 
end 
% transf. of z1~2 to x1~2 
T12 = chol(cov12,'lower');
x(:,1:2) = (T12*x(:,1:2)')';
for i=1:2
    x(:,i) = x(:,i) + NMN(i);
end 

% Covariance matrix for x3~x4 
cov34 = zeros(2,2); %match index by + 2 
for i = 1:2
    for j = i:2
        cov34(i,j) = cor(i+2,j+2)*SIG(i+2)*SIG(j+2);
    end 
end 
for i = 1:2
    for j=i:2
        cov34(j,i) = cov34(i,j);
    end 
end 
% transf. of z3~4 to x3~4
T34 = chol(cov34,'lower');
x(:,3:4) = (T34*x(:,3:4)')';
for i=1:2
    x(:,i+2) = x(:,i+2) + NMN(i+2);
end 
%% transformation for x5 (weibull distribution) 
k = 3.71377203781000;
lamda = 1.10786387179285;
x(:,5) = (lamda)*(-log(1-normcdf(x(:,5)))).^(1/k);

% transf. of z6 to x6 (gumbel dist.) 
beta6 = 0.185642095531828;
mu6 = 0.89284447439388222364137325872783;
x(:,6) = mu6 + beta6*(-log(-log(normcdf(x(:,6)))));

% transf. of z7 to x7 (gumbel dist.)
beta7 = 0.194924200308419;
mu7 = 0.88748669811357634764020434862424;
x(:,7) = mu7 + beta7*(-log(-log(normcdf(x(:,7))))); 
cova = zeros(nd,nd);
cova(1:2,1:2) = cov12;
cova(3:4,3:4) = cov34;
% monomial basis (MNB) 
MNB = zeros(nSample, nA); 
for i=1: nA
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
MNB1 = MNB;
MNB2 = MNB;
MNB1(:,[INDEX1]) = [];
MNB2(:,[INDEX2]) = [];
INFM1 = (ORN1*MNB1')';
INFM2 = (ORN2*MNB2')';

% generate input-output data 
% LY: # of constraints 
rsvl = zeros(nSample, 2);
for L=1:nSample 
        Y1 = responY1(x(L,:), MN);
		Y2 = responY2(x(L,:), MN);
        rsvl(L,1) = Y1;
		rsvl(L,2) = Y2; 
end 
rssc = zeros(nSample, 2, nd, 2);
for r = 1:2 
    for j = 1:2
        for L = 1:nSample
		if (j==1), opt = 'sub1'; INDX = [1,3,4]; end 
		if (j==2), opt = 'sub2'; INDX = [2,3,4]; end 
		scoredk = zeros(nd,1);
        scoredk(INDX,:) = scfcn_gauss(x(L,:), NMN, cova, opt)';
        if (r==1), rssc(L,j,:,1) = rsvl(L,j)*scoredk; end 
        if (r==2), rssc(L,j,:,2) = (rsvl(L,j)^2)*scoredk; end 
        end 
    end 
end 

CFNs = zeros(nA1,2,nd,2);
CFNy = zeros(nA1, 2);
% estimate GPCE coefficients with the standard least sequares
CFNy(:,1) = (INFM1'*INFM1)\(INFM1'*rsvl(:,1));
CFNy(:,2) = (INFM2'*INFM2)\(INFM2'*rsvl(:,2));	
 
for r=1:2
    for j=1:2 
        for i=1:nd 
		if (j==1), CFNs(:,j,i,r) = (INFM1'*INFM1)\(INFM1'*rssc(:,j,i,r)); end 
		if (j==2), CFNs(:,j,i,r) = (INFM2'*INFM2)\(INFM2'*rssc(:,j,i,r)); end  
        end
    end 
end 

%% Sensitivity analysis 
% first moment
jcb = zeros(nd,nd);
for i = 1:nd
    jcb(i,i) = 1/dv(i);
end 

sen1m = zeros(2,nd);
for k=1:2
	for i = 1:nd
		sen1m(k,i) = sen1m(k,i) + CFNs(1,k,i,1);
	end
end 
 
sen1m = sen1m*jcb;
% second moment 
sen2m = zeros(2,nd);
for p = 1:2
    for i = 1:nd 
   sen2m(p,i) = sen2m(p,i) + CFNs(1,p,i,2);
    end 
end 

sen2m = sen2m*jcb;

varaOfcnst  = sum(CFNy(2:end,:).^2);
meaOfcnst = CFNy(1,:);

cnst = (meaOfcnst + 3.*sqrt(varaOfcnst));
cnst1 = cnst(:);
ceq = [];
for i=1:2
	grdnV(:,i) = sen1m(i,:)' + (3/2)*(1/sqrt(varaOfcnst(i)))*(sen2m(i,:)-2*meaOfcnst(i)*sen1m(i,:))';
end 
diffcnst = grdnV(:);
grdneq = [];
disp('Inequality:');
disp(cnst);
save(FilNam2,'x','INFM1','INFM2', 'cova');


else % 2nd or higher iternation 

load(FilNam1); %% load ID, QQ 
load(FilNam2);
% if basis order is not enough, then increase the sampling size than five times of it. 
nSample = nA*15;

% generate input-output data 
rsvl = zeros(nSample, 2);
for L=1:nSample 
        Y1 = responY1(x(L,:), MN);
		Y2 = responY2(x(L,:), MN);
        rsvl(L,1) = Y1;
		rsvl(L,2) = Y2; 
end 

rssc = zeros(nSample, 2, nd, 2);
for r = 1:2 
    for j = 1:2
        for L = 1:nSample
		if (j==1), opt = 'sub1'; INDX = [1,3,4]; end 
		if (j==2), opt = 'sub2'; INDX = [2,3,4]; end 
		scoredk = zeros(nd,1);
        scoredk(INDX,:) = scfcn_gauss(x(L,:), NMN, cova, opt)';
        if (r==1), rssc(L,j,:,1) = rsvl(L,j)*scoredk; end 
        if (r==2), rssc(L,j,:,2) = (rsvl(L,j)^2)*scoredk; end 
        end 
    end 
end 

CFNs = zeros(nA1,2,nd,2);
CFNy = zeros(nA1, 2);
% estimate GPCE coefficients with the standard least sequares 

CFNy(:,1) = (INFM1'*INFM1)\(INFM1'*rsvl(:,1));
CFNy(:,2) = (INFM2'*INFM2)\(INFM2'*rsvl(:,2));	
 
for r=1:2
    for j=1:2 
        for i=1:nd 
		if (j==1), CFNs(:,j,i,r) = (INFM1'*INFM1)\(INFM1'*rssc(:,j,i,r)); end 
		if (j==2), CFNs(:,j,i,r) = (INFM2'*INFM2)\(INFM2'*rssc(:,j,i,r)); end  
        end
    end 
end 
%% Sensitivity analysis 
% first moment
jcb = zeros(nd,nd);
for i = 1:nd
    jcb(i,i) = 1/dv(i);
end 

sen1m = zeros(2,nd);
for k=1:2
	for i = 1:nd
		sen1m(k,i) = sen1m(k,i) + CFNs(1,k,i,1);
	end
end 
 
sen1m = sen1m*jcb;
% second moment 
sen2m = zeros(2,nd);
for p = 1:2
    for i = 1:nd 
   sen2m(p,i) = sen2m(p,i) + CFNs(1,p,i,2);
    end 
end 

sen2m = sen2m*jcb;

varaOfcnst  = sum(CFNy(2:end,:).^2);
meaOfcnst = CFNy(1,:);

cnst = (meaOfcnst + 3.*sqrt(varaOfcnst));

ceq = [];
for i=1:2
	grdnV(:,i) = sen1m(i,:)' + (3/2)*(1/sqrt(varaOfcnst(i)))*(sen2m(i,:)-2*meaOfcnst(i)*sen1m(i,:))';
end 

grdneq = [];
disp('Inequality:');
disp(cnst);

end % end of constraint func. 










            

















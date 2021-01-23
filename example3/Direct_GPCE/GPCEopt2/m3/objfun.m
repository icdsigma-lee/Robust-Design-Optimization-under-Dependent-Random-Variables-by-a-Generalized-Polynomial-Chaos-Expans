%% ========================================================================
% Example 3: Function of objective function and its grandient (Direct GPCE w/partitioned DMORPH)   
% Input: design variables (dv) 
% Output: objective value and design sensitivities
% written by Dongjin Lee (dongjin-lee@uiowa.edu)
%% ========================================================================
function [valObj , vObjGrad] = objfun(dv)

global cntObj stat0 statf sopt 

double precision;
%% Initialization
% number of variables
cntObj = cntObj + 1;
N = 7; % number of random variables 
nd = 4; % # of design variables 
m =3; % order of GPCE for generic function 

% # of GPCE coefficients 
nA = nchoosek(N+m, m);
nA0 = nchoosek(N-2+m,m);

FilNam1 = sprintf('frtst.mat');
FilNam2 = sprintf('secstobj.mat');
% Define parameters of distributiosn
% normalized mean vector for N variables 
NMN = ones(1,N);
tmu5 = 10000; %kgm^3
tmu6 = 2050000000; %Pa
tmu7 = 200000; %N
% transform to original mean vector 
MN = [dv(1), dv(2), dv(3), dv(4), tmu5, tmu6, tmu7];
% standard deviatiaon from 
% SIG(11) is unnecessary 
SIG = [0.02, 0.02, 0.02, 0.02, 0.3, 25/105, 0.25];

% Correlation coefficient matrix in normal dist. 
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
if (cntObj == 1)
load(FilNam1); %% load ID, ORN 

% If basis order is fairly low than that of objective function, then increase sampling size five times more than that size. 
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
cova = zeros(nd,nd);
cova(1:2,1:2) = cov12;
cova(3:4,3:4) = cov34;
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

% Monomial basis (MNB)
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

MNB0 = MNB;
MNB0(:,[INDEX0]) = [];
INFM0 = (ORN0*MNB0')';

% output data   
rsvl = zeros(nSample,1);

% first two order response-values times score function   
rssc = zeros(nSample,nd, 2);  

for L = 1:nSample  
        w = responY0(x(L,:), MN);
        rsvl(L,1) = w;
end 
% 
for r = 1:2 
        for L = 1:nSample
        scoredk = scfcn_gauss(x(L,:), NMN, cova, 'full')';
        if (r==1), rssc(L,:,1) = rsvl(L,1)*scoredk; end % integrand of first moment
        if (r==2), rssc(L,:,2) = (rsvl(L,1)^2)*scoredk; end % integrand of second moment 
        end 
end 
CFNs = zeros(nA0,nd,2);            
CFNy = (INFM0'*INFM0)\(INFM0'*rsvl);
for r = 1:2
    for i=1:nd 
        CFNs(:,i,r) = (INFM0'*INFM0)\(INFM0'*rssc(:,i,r));
    end
end 

%% Sensitivity analysis
% first moment
jcb = zeros(nd,nd);
for i = 1:nd
    jcb(i,i) = 1/dv(i);
end 

sen1m = zeros(1,nd);
for i = 1:nd
    sen1m(1,i) = sen1m(1,i) + CFNs(1,i,1);
end 

sen1m = sen1m*jcb;
% second moment 
sen2m = zeros(1,nd);

for i=1:nd %nA=m
    sen2m(1,i) = sen2m(1,i) + CFNs(1,i,2);
end 
    
sen2m = sen2m*jcb;

varaOfobj = sum(CFNy(2:end).^2);
meaOfobj = CFNy(1);

w1 = meaOfobj;
w2 = sqrt(varaOfobj);
valObj = 0.5*meaOfobj/w1 + 0.5*sqrt(varaOfobj)/w2;
vObjGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varaOfobj))*(sen2m - 2*meaOfobj*sen1m))'; 

disp(dv) 
disp('var of obj:');
disp(varaOfobj)
disp('mean of obj:');
disp(meaOfobj)

% data saving for 2nd step 
save(FilNam2, 'INFM0', 'x', 'w1','w2','cova');

% from 2nd iterations 
else %(cntObj~=1)
    load(FilNam1);
    load(FilNam2); %call INFM, x, w1, w2
    nSample = nA*15;
    % (Part needed for updating at next steps) 
rsvl = zeros(nSample,1);
% first two order response-values times score function   
rssc = zeros(nSample,nd, 2);  
for L = 1:nSample 
     w = responY0(x(L,:), MN);
     rsvl(L,1) = w;
end 

for r = 1:2 
        for L = 1:nSample
        scoredk = scfcn_gauss(x(L,:), NMN, cova, 'full')';
        if (r==1), rssc(L,:,1) = rsvl(L,1)*scoredk; end 
        if (r==2), rssc(L,:,2) = (rsvl(L,1)^2)*scoredk; end 
        end 
end 

CFNs = zeros(nA0,nd,2);
CFNy = (INFM0'*INFM0)\(INFM0'*rsvl);
for r=1:2
    for i=1:nd 
        CFNs(:,i,r) = (INFM0'*INFM0)\(INFM0'*rssc(:,i,r));
    end 
end

%% Sensitivity analysis 
% first moment
jcb = zeros(nd,nd);
for i = 1:nd
    jcb(i,i) = 1/dv(i);
end 

sen1m = zeros(1,nd);
for i = 1:nd
    sen1m(1,i) = sen1m(1,i) + CFNs(1,i,1);
end 

sen1m = sen1m*jcb;

% second moment 
sen2m = zeros(1,nd);
for i=1:nd %nA=m
    sen2m(1,i) = sen2m(1,i) + CFNs(1,i,2);
end    

sen2m = sen2m*jcb;
    
varaOfobj = sum(CFNy(2:end).^2);
meaOfobj = CFNy(1);


valObj = 0.5*meaOfobj/w1 + 0.5*sqrt(varaOfobj)/w2;
vObjGrad = ((0.5/w1)*sen1m + (0.5/(2*w2))*(1/sqrt(varaOfobj))*(sen2m - 2*meaOfobj*sen1m))';  
disp(dv);
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



            

















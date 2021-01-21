%% version 2 : 10 bar truss FEM analysis 
% output : Stresses and displacements 
% written by D. Lee (dongjin-lee@uiowa.edu)

function [YY, LY, w]=RSPSF_TRUSSf2(X, mu_orgn, opt)
global cntRspObj cntRspCon


NX = length(X);

% initialize Z 
Z = zeros(NX,1);

for i = 1:NX
    X(i) = X(i)*mu_orgn(i);
    Z(i) = X(i);
end 

% US unit (in, lb, lbf)  
% Check if input variables are cell or vector 
crt=iscell(X); %0: Vector, 1: Cell 
w = 0; % initialize weight of truss 
switch crt
%% CASE: VECTOR     
    case 0           
% information of truss 
nElem=10; 
nDof=2; % d.o.f of bar element 
nNode=6; % # of nodes 
PV=-1E5; % Unit: lbf, vertical loading 
PH =4E5; % Unit: lbf, horizontal loading 
L=360; % Unit: in 
rho = 0.1; % Unit: lb/in^3 
E=10^7; % Unit: lbf/in^2  
C=E/L; 

% initialize truss FEM  
elem=zeros(nElem,2); % element matrix 
KK=zeros(nDof*nNode,nDof*nNode); % stiffness matrix
FF=zeros(nDof*nNode,1); % external force vector
UU=zeros(nDof*nNode,1); % displancement vector 
% elem = [node1 node2] 
elem(1,:)=[1 2];
elem(2,:)=[2 3];
elem(3,:)=[3 6];
elem(4,:)=[5 6];
elem(5,:)=[4 5];
elem(6,:)=[1 5];
elem(7,:)=[2 4];
elem(8,:)=[2 5];
elem(9,:)=[2 6];
elem(10,:)=[3 5];
a1=cos(pi/4); a2=sin(pi/4);
b1=-cos(pi/4); b2=sin(pi/4);
% Transformation matrix 
A=[a1 a2 0 0;-a2 a1 0 0;0 0 a1 a2;0 0 -a2 a1];
B=[b1 b2 0 0;-b2 b1 0 0;0 0 b1 b2;0 0 -b2 b1];

% Generating Global Stiffness Matrix 
switch opt
    case 'obj'
cntRspObj = cntRspObj + 1;       
for ii=1:nElem 
    switch ii 
        case {1, 2, 4, 5}     
             w = w + Z(ii)*L;
        case {3, 8}
             w = w + Z(ii)*L;
        case {6, 9} 
             w = w + Z(ii)*L*sqrt(2);
        case {7, 10}
             w = w + Z(ii)*L*sqrt(2);
        otherwise 
            disp('wrong element');
    end 
end

% weight 
w = rho*w; % Unit: lb 
YY = 0; % default of perform. func. value 
LY = 0; % # of perform. funct. 
    case 'cnstn' %constraint function
        cntRspCon = cntRspCon + 1;
%% Construct global stiffness matrix 
        for ii=1:nElem 
    index=[elem(ii,1)*nDof-1:elem(ii,1)*nDof, elem(ii,2)*nDof-1:elem(ii,2)*nDof];
    switch ii 
        case {1, 2, 4, 5}
            % Z(ii) : input random values of area 
            KK(index,index)=KK(index,index)+Z(ii)*C*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0];
        case {3, 8}
            KK(index,index)=KK(index,index)+Z(ii)*C*[0 0 0 0;0 1 0 -1;0 0 0 0;0 -1 0 1];
        case {6, 9}
            KK(index,index)=KK(index,index)+Z(ii)*(C/sqrt(2))*A'*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0]*A;
        case {7, 10}
            KK(index,index)=KK(index,index)+Z(ii)*(C/sqrt(2))*B'*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0]*B;
        otherwise 
            disp('wrong element');
    end 
    end 
    
% Generation of Force matrix 
FF(2*nDof) = PV;
FF(3*nDof) = PV;
FF(3*nDof-1)= PH;
% Boundary condtion 
    KK([1*nDof-1:1*nDof,4*nDof-1:4*nDof], :)=[];
    KK(:,[1*nDof-1:1*nDof,4*nDof-1:4*nDof])=[];
    FF([1*nDof-1:1*nDof,4*nDof-1:4*nDof])=[];
% Solve displacement  
U=inv(KK)*FF;

sIndex=[2 3 5 6];
nS=length(sIndex);
for jj=1:nS 
    UU(sIndex(jj)*nDof-1:sIndex(jj)*nDof)=U(jj*nDof-1:jj*nDof);
end 
% end of displacement solution 
%% Construct stress 
con = sqrt(2)/2;
S = zeros(nElem, 1);
for ii=1:nElem 
    index=[elem(ii,1)*nDof-1:elem(ii,1)*nDof,elem(ii,2)*nDof-1:elem(ii,2)*nDof];
    switch ii 
        case {1, 2, 4, 5}
            S(ii,1)=C*[-1 0 1 0]*UU(index);
        case {3, 8}
            S(ii,1)=C*[0 -1 0 1]*UU(index); 
        case {6, 9}
            S(ii,1)=(C/sqrt(2))*[-con -con con con]*UU(index); 
        case {7, 10}
            S(ii,1)=(C/sqrt(2))*[con -con -con con]*UU(index);
        otherwise 
            disp('wrong element');
    end 
end 
% define performance function 
LU = 1;
LS = nElem;
LY = LU + LS;
YY = zeros(1, LY);
YY(1, 1:LU) = (-1 + abs(UU(3*nDof))/5)'; % unit, in
for i=1:nElem
	    %YY(1,i+1) = S(i);
		if (i==10)
        YY(1,i+1) = (-1 + abs(S(i))/75000); %if (i=10)
        else 
		YY(1,i+1) = (-1 + abs(S(i))/25000);  
		end 
end 

    otherwise 
    disp('wrong input');
end 

otherwise end 
        

    

              
            

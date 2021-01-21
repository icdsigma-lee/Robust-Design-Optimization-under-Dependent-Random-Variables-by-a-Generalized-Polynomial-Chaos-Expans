%% ========================================================================
%  submain program (example4) Multi-Point Single-Step scheme 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [history,searchdir, subUb, subLb, df] = runfmincon(x0)
diary historywindow
diary on 
global cntObj cntCon beta ub lb ii
filNam = sprintf('sizing.mat');
history.x = [];
history.fval = [];

searchdir = [];
nd = 10;
fir = 0; % initial iteration 
% golden ratio constant 
gr = (1+sqrt(5))/2; 
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','OutputFcn',@outfun);
options = optimoptions(options,'GradObj','on','GradConstr','on','TolCon',1e-4,'TolFun',1e-4,'TolX',1e-4);
lb = ones(1, nd)*0.1; ub = ones(1,nd)*35; 
cntObj = 1; cntCon = 1;  

% terminiation condition
eps1 = 1E-6; % stopping criteria for design variables
eps3 = 1E-6; % stopping criteria for objective value 

% feasible estimation 
dfLast = (ub - x0);
q = 0;
objValue = 1;
objValueLast = 1+1; % 

% sizing factor of subdomain (default: 0.3) 
beta = ones(nd,1)*0.3;

cnt = 0; % if tcnt == 0, Do the below loop, and otherwise, stop! 
objValueLast0 = 10;
df = zeros(1,nd);
while (cnt == 0)
% q iteration
q = q + 1;

cntSizing = 0; %
cntSizing1 = 0;
tmp = 0; 
while (cntSizing == 0)
	tmp = tmp+1;
	cntCon = 0; cntObj = 0; ii=0;
    if (fir ==0), cntCon = 1; cntObj = 1; end 
	[objValue, ~] = sobjfun(x0);
	[c, ~, ~, ~] = sconfun(x0);
    fir = fir + 1;
	
	chk = 0;	
	for i=1:length(c)
		if (c(i) > 1e-1), chk = 1; end 
	end 
	if (chk == 1)
	cntSizing1 = cntSizing1 + 1;
	x0 = dfLast*(1/gr) + x0*(1-1/gr); % find feasible design
	df(q,:) = x0;
	else 
	df(q,:) = x0;
	cntSizing = 1;
	end 
	if (cntSizing1 == 10)
		disp('Searching feasible design fails');
		cntSizing = 1;
	end 
	
% sizing of subdomain
if ((tmp==1) && (q~=1)) 
	err1 = sqrt((objValue - objValueLast)'*(objValue - objValueLast));
    err2 = sqrt((c-cLast)*(c-cLast)');
	if ((err1 <= 1e-2) && (err2 <= 1e-2))
	% increase overall subregion size 
	beta(:) = beta(:).*gr;
	elseif ((err1 >= 7e-2) || (err2 >= 7e-2))
	% decrease overall subregion size 
	beta(:) = beta(:).*1/gr;
    else end 
	for i=1:nd
		% Is constraint inequality active?
		chkL = sqrt((x(i)-subLb(q-1,i))^2);
		chkU = sqrt((x(i)-subUb(q-1,i))^2);
		if ((chkL < 1e-2) || (chkU < 1e-2))
            beta(i) = beta(i).*gr;
		end 
		% Is dv on convergence step 
        err3 = sqrt((df(q,i)-dfLast(i))^2);
		if (err3 < 0.5)
		beta(i) = beta(i).*(1/gr);
        end 
        % limit subregion size
		if (beta(i) < 0.05)
			beta(i) = 0.05;
		end 
    end
end 
end 

% decide subregion
for i=1:nd 
    %active check 
    if (q ~= 1)
     if (abs(x(i)-subLb(q-1,i)) < 1e-4)    
     subLb(q,i) = df(q,i)  -  beta(i)*(ub(i) - lb(i))*gr/(gr+1);
     subUb(q,i) = df(q,i)  +  beta(i)*(ub(i) - lb(i))/(gr+1);
     else 
     subLb(q,i) = df(q,i) -   beta(i)*(ub(i) - lb(i))/2;
     subUb(q,i) = df(q,i)  +  beta(i)*(ub(i) - lb(i))/2;
     end 
     
     if (abs(x(i)-subUb(q-1,i)) < 1e-4)    
     subLb(q,i) = df(q,i)  -  beta(i)*(ub(i) - lb(i))/(gr+1);
     subUb(q,i) = df(q,i)  +  beta(i)*(ub(i) - lb(i))*gr/(gr+1);
     else 
     subLb(q,i) = df(q,i) -   beta(i)*(ub(i) - lb(i))/2;
     subUb(q,i) = df(q,i)  +  beta(i)*(ub(i) - lb(i))/2;
     end
    else 
        
     subLb(q,i) = df(q,i) -   beta(i)*(ub(i) - lb(i))/2;
     subUb(q,i) = df(q,i)  +  beta(i)*(ub(i) - lb(i))/2;
    end     
     
    if (subLb(q,i) < lb(i))
        subLb(q,i) = lb(i);
    end 
    if (subUb(q,i) > ub(i))
        subUb(q,i) = ub(i);
    end         
end

disp('sub-upper');
disp(subUb(q,:));
disp('sub-lower');
disp(subLb(q,:));	
save(filNam, 'subUb', 'subLb');
terminate1 = (df(q,:)-dfLast)*(df(q,:)-dfLast)';
terminate3 = abs((objValue - objValueLast0) / objValue); 
objValueLast0 = objValue;	

cnt = 0; 
if (terminate1 < eps1), cnt = 1; end 
if (terminate3 < eps3), cnt = 1; end 

switch (cnt)
    case 1
        disp('Optimization is completed');
    case 0      
        cntCon = 1; cntObj= 1; % single-step 
        % solving local RDO problem using single-step GPCE
        [x,fval] = fmincon(@sobjfun,df(q,:),[],[],[],[],subLb(q,:),subUb(q,:),@sconfun,options);
        % obtaining data from step2 
        [objValueLast,~] = sobjfun(x);
        [cLast, ~, DCLast, ~] = sconfun(x);
    otherwise  
         disp('Option was wrongly selected');
end 
dfLast = df(q,:);
x0 = x;
end 

% define outfun to save information at each iteration 
function stop = outfun(x,optimValues,state)
stop = false;
FilNam = sprintf('backup.mat');
 
   switch state
       case 'init'
           hold on
       case 'iter'
           % Concatenate current point and objective function
           % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
           % Concatenate current search direction with 
           % searchdir.
           searchdir = [searchdir;...
                        optimValues.searchdirection'];
		   save(FilNam,'history');
%            plot(x(1),x(2),'o');
%            % Label points with iteration number.
%            % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,x(2),num2str(optimValues.iteration));
       case 'done'
           hold off
       otherwise
   end
end
diary off 
end 

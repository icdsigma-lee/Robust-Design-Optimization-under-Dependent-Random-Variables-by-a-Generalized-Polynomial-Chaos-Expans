%% ========================================================================
%  submain program (example4)
%  run optimization 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [history,searchdir] = runfmincon(x0)
global cntObj cntCon cntRspObj cntRspCon
history.x = [];
history.fval = [];
searchdir = [];
nd = 10;

options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','OutputFcn',@outfun);
options = optimoptions(options,'GradObj','on','GradConstr','on','TolCon',1e-5,'TolFun',1e-4,'TolX',1e-5);
lb = ones(1, nd)*0.1; % lower bounds 
ub = ones(1,nd)*35;   % upper bounds
cntObj = 1; cntCon = 1; cntRspObj = 0; cntRspCon = 0;

[x,fval] = fmincon(@objfun,x0,[],[],[],[],lb,ub,... 
   @confun,options);

% define outfun to save information at each iteration 
function stop = outfun(x,optimValues,state)
stop = false;
 
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
       case 'done'
           hold off
       otherwise
   end
end

end 
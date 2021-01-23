%% ========================================================================
%  submain program (example1-2)
%  run optimization 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================

function [history,searchdir] = runfmincon(x0)
history.x = [];
history.fval = [];
searchdir = [];
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','OutputFcn',@outfun);
options = optimoptions(options,'GradObj','on','GradConstr','on','TolCon',1e-10,'TolFun',1e-10,'TolX',1e-15);
lb = [0,0];     % lower bounds 
ub = [10,10];   % upper bounds

% SQP optimization 
% x0 : initial design
[x,fval] = fmincon(@exactObjfun,x0,[],[],[],[],lb,ub,... 
   @exactConfun,options);

% save information at every iteration  
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
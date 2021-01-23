%% ========================================================================
%  submain program (example1-1)
%  run optimization 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================

function [history,searchdir] = runfmincon(x0)
global cntObj cntCon cntY1 cntY2 
history.x = [];
history.fval = [];
searchdir = [];
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','OutputFcn',@outfun);
options = optimoptions(options,'GradObj','on','GradConstr','on');
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
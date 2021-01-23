%% ========================================================================
%  submain program (example3)
%  run optimization 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [history,searchdir] = runfmincon(x0)
global cntObj cntCon con  cntY0 cntY1 cntY2 tGram
history.x = [];
history.fval = [];
searchdir = [];
con = zeros(2,1);
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','OutputFcn',@outfun);
options = optimoptions(options,'GradObj','on','GradConstr','on');
lb = [2,2,0.3,0.3]; % lower bounds
ub = [25,25,1.4,1.4];   % upper bounds
cntObj = 1; cntCon = 0; cntY0 = 0; cntY1 = 0; cntY2 = 0; tGram = 0;

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
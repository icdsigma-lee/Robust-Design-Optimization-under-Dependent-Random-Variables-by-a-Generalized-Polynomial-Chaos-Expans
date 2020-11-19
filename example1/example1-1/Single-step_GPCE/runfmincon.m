%% ========================================================================
%  submain program (example1-1)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [history,searchdir] = runfmincon(x0)
global cntObj cntCon cntY1 cntY2 
history.x = [];
history.fval = [];
searchdir = [];
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','OutputFcn',@outfun);
options = optimoptions(options,'GradObj','on','GradConstr','on');
lb = [0,0];  ub = [10,10];   % No upper or lower bounds

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
%            plot(x(1),x(2),'o');
%            % Label points with iteration number.
%            % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,x(2),num2str(optimValues.iteration));
       case 'done'
           hold off
       otherwise
   end
end

end 
%% ========================================================================
%  score function (Multivariate Gaussian) 
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%% ========================================================================
function [output]=scoreGauss(X, mu, cov, opt)
%Multivariate Gaussian score function  
mu1 = mu(1);
mu2 = mu(2);
mu3 = mu(3);
mu4 = mu(4);
switch opt
    case 'full' 
% X1-X4
dist = [mu1*X(1)-mu1,mu2*X(2)-mu2,mu3*X(3)-mu3,mu4*X(4)-mu4]; 
output = cov\dist';
    case 'sub1'
        
dist = [mu1*X(1)-mu1,mu3*X(3)-mu3,mu4*X(4)-mu4];
cov(:,2) = [];
cov(2,:) = [];
output = cov\dist';

    case 'sub2'
        
dist = [mu2*X(2)-mu2,mu3*X(3)-mu3,mu4*X(4)-mu4];
cov(:,1) = [];
cov(1,:) = [];
output = cov\dist';
    otherwise 
end 


end 
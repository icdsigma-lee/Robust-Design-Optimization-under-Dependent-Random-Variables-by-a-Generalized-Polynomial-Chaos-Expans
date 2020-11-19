% Function "bas_h" outputs monomials for Hermite polynomials {x^(j)} 
function [outPut] = bas(x,j)
nTemp = length(x);
% Find mean value when j is larger than 0 and even number 
  if (nTemp == 1)
     outPut = x^(j);
 elseif (nTemp > 1)
     outPut = x.^(j);
 end 

end 


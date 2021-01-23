function [outPut] = bas(x,j)
nTemp = length(x);
  if (nTemp == 1)
     outPut = x^(j);
 elseif (nTemp > 1)
     outPut = x.^(j);
 end 

end 


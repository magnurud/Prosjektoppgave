function [Ah] = stiffness_2D(N,x,y,wX,wY,LDM);
% description:
%      generate the stiffness matrix Ah;
%
% arguments:
%   - N     Number of points in each direction  
%   - x     x values of the coordinates 
%   - y     y values of the coordinates 
%   - wX    The weight used in the quadrature rule for the x-direction
%   - wY    The weight used in the quadrature rule for the x-direction
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
% returns:
%		- Ah 	  Stiffness matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: April 2015

dofs = N^2;
Ah = zeros(dofs);

for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  for J = 1:dofs
    k = mod(J-1,N)+1;
    l = fix((J-1)/N)+1;
    if(l==j) 
      for alpha = 1:N
        %Ah(I,J) = Ah(I,J) + wX(alpha)*wY(l)*LDM(i,alpha)*LDM(k,alpha);
        Ah(I,J) = Ah(I,J) + wX(alpha)*wY(l)*LDM(alpha,i)*LDM(alpha,k);
      end
    end
    if(i==k) 
      for beta = 1:N
        %Ah(I,J) = Ah(I,J) + wX(i)*wY(beta)*LDM(j,beta)*LDM(l,beta);
        Ah(I,J) = Ah(I,J) + wX(i)*wY(beta)*LDM(beta,j)*LDM(beta,l);
      end
    end
  end
end

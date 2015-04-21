function [Gh] = gradient_2D(LDM,B1,B2,N,wX,wY);
% description:
%      generate the gradient matrix Gh;
%
% arguments:
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
%   - B1     The N by N matrix containing the first component of the vectorfield B evaluated in each
%   (x,y)
%   - B2     The N by N matrix containing the first component of the vectorfield B evaluated in each
%   (x,y)
%   - N     Number of points in each direction  
%   - wX    The weight used in the quadrature rule for the x-direction
%   - wY    The weight used in the quadrature rule for the y-direction
% returns:
%		- Ah 	  Stiffness matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: April 2015

dofs = N^2;
Gh = zeros(dofs);
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  for J = 1:dofs
    k = mod(J-1,N)+1;
    l = fix((J-1)/N)+1;
    if(l==j) 
        Gh(I,J) = Gh(I,J) + wX(i)*wY(j)*B1(i,j)*LDM(i,k);
    end
    if(i==k) 
        Gh(I,J) = Gh(I,J) + wX(i)*wY(j)*B2(i,j)*LDM(j,l);
    end
  end
end


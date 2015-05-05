function [Jfh] = load_LS(N,x,y,wX,wY,f,LDM,dB1,dB2)
%function [Jfh] = load_2D_fast(N,x,y,wX,wY,f,LDM)
%
% description:
%     calculates the inner product <f,v> on the unit square using Spectral
%     LS methods.
%
% arguments:
%   - N     Number of points in each direction. 
%   - x     x-values of the coordinates
%   - y     y-values of the coordinates
%   - wX    quadrature weights
%   - wY    quadrature weights 
%		- f 		Loading function expression of two variables
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
%		- dB1		The derivative of the first component of the vector field wrt u diagonal matrix
%		- dB2		The derivative of the second component of the vector field wrt u
% returns:
%   - Jfh 		The jacobian of the loading function evaluated in all the nodes.
%
% author: Magnus Aa. Rud
% last edit: April 2015

NLS = 3*N; 
dofs = N^2;
LSdofs = 3*dofs;
Jfh = zeros(LSdofs,1);
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  Jfh(I) = Jfh(I)           - wX(i)*wY(j)*dB1(i,j)*f(x(i),y(j));
  Jfh(I+dofs) = Jfh(I+dofs) - wX(i)*wY(j)*dB2(i,j)*f(x(i),y(j));
end

Jfh = diag(Jfh);

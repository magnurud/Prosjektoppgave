function [fh] = load_2D_precond(N,x,y,wX,wY,f,LDM)
%function [fh] = load_2D_precond(N,x,y,wX,wY,f,LDM)
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
% returns:
%   - fh 		The loading function evaluated in all the nodes.
%
% author: Magnus Aa. Rud
% last edit: April 2015

NLS = 3*N; 
dofs = N^2;
LSdofs = 3*dofs;
fh = zeros(LSdofs,1);
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
	for gamma = 1:N
  fh(I) = fh(I) + wX(gamma)*wY(j)*LDM(gamma,i)*f(x(gamma),y(j));
  fh(I+dofs) = fh(I+dofs) + wX(i)*wY(gamma)*LDM(gamma,j)*f(x(i),y(gamma));
  end
end



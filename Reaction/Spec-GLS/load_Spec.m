function [fh] = load_2D(N,x,y,wX,wY,f)
%function [fh] = load_2D(N,x,y,wX,wY,f)
%
% description:
%     calculates the inner product <f,v> on the unit square using Spectral
%     methods.
%
% arguments:
%   - N     Number of points in each direction. 
%   - x     x-values of the coordinates
%   - y     y-values of the coordinates
%   - wX    quadrature weights
%   - wY    quadrature weights 
%		- f 		Loading function expression of two variables
% returns:
%   - fh 		The loading function evaluated in all the nodes.
%
% author: Magnus Aa. Rud
% last edit: April 2015

fh = zeros(N^2,1);
for I = 1:N^2
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  fh(I) = wX(i)*wY(j)*f(x(i),y(j));
end


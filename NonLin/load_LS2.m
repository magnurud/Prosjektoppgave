function [fh] = load_LS2(W,LDM,B1,B2,F,mu,sigma)
%function [fh] = load_2D_fast(N,x,y,wX,wY,f,LDM)
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
%   - mu    The diffusion constant
%   - sigma The reaction constant
% returns:
%   - fh 		The loading function evaluated in all the nodes.
%
% author: Magnus Aa. Rud
% last edit: April 2015

dofs = length(F);
fh = zeros(3*dofs,1);

fh(1:dofs) = mu*(kron(W,LDM'*W))*F - (kron(W,W))*B1*F; 

fh(dofs+1:2*dofs) = mu*(kron(LDM'*W,W))*F-(kron(W,W))*B2*F; 

%fh(2*dofs+1:end) =  sigma*kron(W,W)*F;


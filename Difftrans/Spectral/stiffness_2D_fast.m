function [Ah] = stiffness_2D_fast(wX,LDM);
% description:
%      generate the stiffness matrix Ah;
%
% arguments:
%   - wX    The weight used in the quadrature rule for the x-direction
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
% returns:
%		- Ah 	  Stiffness matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: April 2015

W = diag(wX);
Ah = kron(W,LDM'*W*LDM)+kron(LDM'*W*LDM,W);

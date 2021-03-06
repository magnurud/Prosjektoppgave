function [Ah] = stiffness_Spec(W,LDM);
% description:
%      generate the stiffness matrix Ah;
%
% arguments:
%   - W    The weight used in the quadrature rule 
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
% returns:
%		- Ah 	  Stiffness matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: April 2015

Ah = kron(W,LDM'*W*LDM)+kron(LDM'*W*LDM,W);

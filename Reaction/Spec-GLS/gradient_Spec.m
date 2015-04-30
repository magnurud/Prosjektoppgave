function [Gh] = gradient_Spec(LDM,B1,B2,W,dofs);
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
%   - W     The GLL-weights
%   - dofs  degrees of freedom 
% returns:
%		- Gh 	  gradient matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: April 2015

Bd1 = diag(reshape(B1,dofs,1)); 
Bd2 = diag(reshape(B2,dofs,1)); 

Gh = Bd1*kron(W,W*LDM) + Bd2*kron(W*LDM,W);


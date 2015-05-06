function [Jh] = jacobi_Spec_lin(W,LDM,dB1,dB2,B1,B2,U,Rg);
% description:
%      generate the jacobian for the LS-part of the nonlin system 
%
% arguments:
%   - W     The matrix with the GLL-weights along the diagonal
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
%		- dB1		The derivative of the first component of the vector field wrt u diagonal matrix
%		- dB2		The derivative of the second component of the vector field wrt u
%		- B1		the first component of the vector field diagonal matrix
%		- B2		the second component of the vector field diagonal matrix
%		- U			the current solution, vector dofs x 1
%		- Rg 		the lifting function, vector dofs x 1
% returns:
%		- Jh 	  Jacobi matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: May 2015

%%% OOPS DIAGONAL MATRICES ARE REALLLY SLOW!!!%%% 
PHI = kron(W,W*LDM);
PSI = kron(WL,W);

Jh = B1*PHI+diag(dB1*PHI*(U+Rg))+B2*PSI+diag(dB2*PSI*(U+Rg)); % From the grad-part Includes BC's


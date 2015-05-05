function [Jh] = jacobi_LS_lin(W,LDM,dB1,dB2,B1,B2,U,W1,mu);
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
%		- W1			the current solution of the first flux component, dofs x 1
%   - mu    the diffusion constant
% returns:
%		- Jh 	  Jacobi matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: May 2015

%%% OOPS DIAGONAL MATRICES ARE REALLLY SLOW!!!%%% 
PHI = kron(W,W*LDM);
PSI = kron(WL,W);
WW = kron(WW);

DG11 = -mu*(B1*PHI+diag(dB1*PHI*U)) -mu*(PHI'*B1+diag(dB1*PHI'*U))+diag(2*dB1*WW*B2*U)+B1*WW*B1...
			 -mu*(diag(dB2*PHI*W1)) -mu*diag(dB1*PSI'*W1)+diag(dB1*WW*B2*W1+dB2*WW*B1*W1);

DG21 = -mu*(PHI'*B2+diag(dB2*PHI'*U)) -mu*(B1*PSI+diag(dB1*PHI*U))+diag(dB2*WW*B1*U+dB1*WW*B2*U)+B2*WW*B2...
			 -mu*(diag(dB2*PSI*W1)) -mu*diag(dB2*PSI'*W1)+diag(2*dB2*WW*B2*W1);

DG12 = -mu*B2*PHI-mu*PSI'*B1+B1*WW*B2;
DG22 = -mu*B2*PSI-mu*PSI'*B2+B2*WW*B2;

ZEROS = zeros(size(WW));


Jh = [DG11 DG12 ZEROS; DG21 DG22 ZEROS; ZEROS ZEROS ZEROS];




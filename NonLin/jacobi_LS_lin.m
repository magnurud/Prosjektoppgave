function [Jh] = jacobi_LS_lin(W,LDM,dB1,dB2,B1,B2,W1,W2,mu);
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
PSI = kron(W*LDM,W);
WW = kron(W,W);


DG13 =  -mu*(diag(dB1*PHI*W1)) -mu*diag(dB1*PHI'*W1) + 2*diag(dB1*WW*B1*W1); % FROM G11*W1
DG13 = DG13 -mu*(diag(dB2*PHI*W2)) -mu*diag(dB1*PSI'*W2)+diag(dB1*WW*B2*W2+dB2*WW*B1*W2); % FROM G12*W2

DG23 = -mu*(diag(dB1*PSI*W1)) - mu*diag(dB2*PHI'*W1) + diag(dB2*WW*B1*W1+dB1*WW*B2*W1); %FROM G21*W1
DG23 = DG23 -mu*(diag(dB2*PSI*W2)) -mu*diag(dB2*PSI'*W2)+diag(2*dB2*WW*B2*W2); %FROM G22*W2


%DG11 = -mu*(B1*PHI+diag(dB1*PHI*U)) -mu*(PHI'*B1+diag(dB1*PHI'*U))+diag(2*dB1*WW*B2*U)+B1*WW*B1...
%			 -mu*(diag(dB2*PHI*W1)) -mu*diag(dB1*PSI'*W1)+diag(dB1*WW*B2*W1+dB2*WW*B1*W1);

%DG21 = -mu*(PHI'*B2+diag(dB2*PHI'*U)) -mu*(B1*PSI+diag(dB1*PHI*U))+diag(dB2*WW*B1*U+dB1*WW*B2*U)+B2*WW*B2...
%			 -mu*(diag(dB2*PSI*W1)) -mu*diag(dB2*PSI'*W1)+diag(2*dB2*WW*B2*W1);

%DG12 = -mu*B2*PHI - mu*PSI'*B1 + B1*WW*B2;
%DG22 = -mu*B2*PSI - mu*PSI'*B2 + B2*WW*B2;

ZEROS = zeros(size(WW));

Jh = [ZEROS ZEROS DG13; ZEROS ZEROS DG23; ZEROS ZEROS ZEROS];




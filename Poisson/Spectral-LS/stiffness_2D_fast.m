function [Ah] = stiffness_2D_fast(W,LDM);
% description:
%      generate the stiffness matrix Ah;
%      This will be organized differently!! 
%
% arguments:
%   - W     The matrix with the GLL-weights along the diagonal
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
% returns:
%		- Ah 	  Stiffness matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: April 2015

A11 = kron(W,LDM'*W*LDM)+kron(W,W);

A12 = kron((W*LDM),(W*LDM)');

A13 = kron(W,W*LDM);

A22 = kron(LDM'*W*LDM,W)+kron(W,W);

A23 = kron(W*LDM,W);

A33 = kron(LDM'*W*LDM,W)+kron(W,LDM'*W*LDM);

Ah = [A11 A12 A13; A12' A22 A23; A13' A23' A33];

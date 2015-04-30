function [Rh] = reaction_LS(mu,B1,B2,sigma,W,LDM,dofs);
% description:
%      generate the reaction matrix Rh;
%
% arguments:
%   - mu    the diffusion constant
%   - B1     The N by N matrix containing the first component of the vectorfield B evaluated in each
%   (x,y)
%   - B2     The N by N matrix containing the first component of the vectorfield B evaluated in each
%   (x,y)
%   - sigma the weight of the reaction term
%   - W     The GLL-weights
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
%   - dofs  degrees of freedom 
% returns:
%		- Rh    The matrix for the reaction term
%
% author: Magnus Aa. Rud
% last edit: April 2015

Bd1 = diag(reshape(B1,dofs,1)); 
Bd2 = diag(reshape(B2,dofs,1)); 

Rh13 = mu*kron(W,LDM'*W) -Bd1*kron(W,W);
Rh23 = mu*kron(LDM'*W,W) -Bd2*kron(W,W);
Rh33 = sigma*kron(W,W);
ZERO = zeros(size(Rh13));

Rh = sparse(sigma*[ZERO ZERO Rh13; ZERO ZERO Rh23 ; Rh13' Rh23' Rh33]);

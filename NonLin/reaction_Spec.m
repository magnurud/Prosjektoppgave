function [Rh] = reaction_Spec(W);
% description:
%      generate the reaction matrix Rh;
%
% arguments:
%   - W     The GLL-weights
% returns:
%		- Rh    The matrix for the reaction term
%
% author: Magnus Aa. Rud
% last edit: April 2015

Rh = kron(W,W);



function [fh] = load_2D(W,F)
%function [fh] = load_2D(N,x,y,wX,wY,f)
%
% description:
%     calculates the inner product <f,v> on the unit square using Spectral
%     methods.
%
% arguments:
%   - W 	The diagonal matrix with the GLL-weights along the diag
%   - F 		Loading function with values along the diag
% returns:
%   - fh 		The loading function evaluated in all the nodes.
%
% author: Magnus Aa. Rud
% last edit: April 2015

fh = kron(W,W)*F;




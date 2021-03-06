function [p tri edge] = getSquare(N),
% function [p tri edge] = getSquare(N),
% 
% description:
%      generate a mesh triangulation of the square (0,1)^2
%
% arguments:
%   - N    the number of nodes in each spatial direction (N^2 total nodes)
% returns:
%   - p     nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   elements. Index to the three corners of element i given in row i.
%   - edge  index list of all nodal points on the outer edge (r=1)
%
% author: Magnus Aa. Rud
% last edit: March 2015

p = zeros(N,2);
p(1,:) = [0,0];

k = 1;
for j=1:N,
	for i=1:N,
		p(k,:) = [(i-1)/(N-1), (j-1)/(N-1)];
		k = k+1;
	end
end

tri  = delaunay(p(:,1), p(:,2)); % Creating the elements
edge = [      1:N];             % south edge
edge = [edge, (2*N):N:(N^2-N)]; % east edge
edge = [edge, N^2:-1:(N^2-N+2)];   % north edge
edge = [edge, (N^2-N+1):-N:(N+1)]; % west edge


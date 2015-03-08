function [p tetr edge] = getBeef(n1,n2,n3),
% function [p tetr edge] = getBeef(n1,n2,n3),
% 
% description:
%      generate a mesh triangulation of a beef. It is described as a parametric
%      volume with the first parametric direction being the length (~15cm),
%      the second being the width (~8cm) and the third being the thickness
%      (~3cm)
%
% arguments
%   - n1   the number of nodes in the first parametric direction
%   - n2   the number of nodes in the second parametric direction
%   - n3   the number of nodes in the third parametric direction
% returns:
%   - p     nodal points. (x,y,z)-coordinates for point i given in row i.
%   - tetr  elements. Index to the four corners of element i given in row i.
%   - edge  index list of all nodal points on the outer edge. Each edge node
%           is also identified with a tag telling which side it originates 
%           from
%             1 :   xi1 = min
%             2 :   xi1 = max
%             3 :   xi2 = min
%             4 :   xi2 = max
%             5 :   xi3 = min
%             6 :   xi3 = max

% author: Kjetil A. Johannessen
% last edit: October 2010

load beef_data
n_basis   = zeros(3,1);
n_basis(1)= length(knot.xi)   - p(1) - 1;
n_basis(2)= length(knot.eta)  - p(2) - 1;
n_basis(3)= length(knot.zeta) - p(3) - 1;

% scaling to get the right size - longest side ~0.15m
scale = eye(3)/40*15/100; 
B = B*scale;

% typo in the problem text, where first and second parametric direction was swapped. Fixing it here instead of in the problem set
tmp = n1;
n1 = n2;
n2 = tmp;

dxi   = knot.xi(end)-knot.xi(1);
deta  = knot.eta(end)-knot.eta(1);
dzeta = knot.zeta(end)-knot.zeta(1);

xi   = knot.xi(1):1/(n1-1)*dxi:knot.xi(end);
eta  = knot.eta(1):1/(n2-1)*deta:knot.eta(end);
zeta = knot.zeta(1):1/(n3-1)*dzeta:knot.zeta(end);

N_xi   = getBSplineBasisAndDerivative(p(1), xi, knot.xi);
N_eta  = getBSplineBasisAndDerivative(p(2), eta, knot.eta);
N_zeta = getBSplineBasisAndDerivative(p(3), zeta, knot.zeta);

p = zeros(n1*n2*n3, 3);

ii = 1;
for zeta=1:n3,
for eta=1:n2,
for xi=1:n1,
	jj = 1;
	for k=1:n_basis(3),
	for j=1:n_basis(2),
	for i=1:n_basis(1),
		p(ii,:) = p(ii,:) + N_xi(i,xi)*N_eta(j,eta)*N_zeta(k,zeta) * B(jj,:);
		jj = jj + 1;
	end
	end
	end
	ii = ii + 1;
end
end
end

tetr = zeros( (n1-1)*(n2-1)*(n3-1)*6, 4);
tetr_i = 1;
for zeta=1:n3-1,
for eta=1:n2-1,
for xi=1:n1-1,
	ii = (zeta-1)*(n1*n2) + (eta-1)*n1 + xi;
	i1 = ii;
	i2 = ii+1;
	i3 = ii+n1;
	i4 = ii+n1+1;
	i5 = ii+n1*n2;
	i6 = ii+n1*n2+1;
	i7 = ii+n1*n2+n1;
	i8 = ii+n1*n2+n1+1;
	tetr(tetr_i, :) = [i1, i2, i3, i7];
	tetr_i = tetr_i + 1;
	tetr(tetr_i, :) = [i2, i5, i6, i7];
	tetr_i = tetr_i + 1;
	tetr(tetr_i, :) = [i1, i2, i5, i7];
	tetr_i = tetr_i + 1;
	tetr(tetr_i, :) = [i2, i6, i7, i8];
	tetr_i = tetr_i + 1;
	tetr(tetr_i, :) = [i2, i3, i4, i7];
	tetr_i = tetr_i + 1;
	tetr(tetr_i, :) = [i2, i4, i7, i8];
	tetr_i = tetr_i + 1;
end
end
end

edge = zeros(n1*n2*2 + n1*n3*2 + n2*n3*2, 2);

edge_i = 1;
glob_i = 1;
for k=1:n3,
for j=1:n2,
for i=1:n1,
	if i==1,
		edge(edge_i,:) = [glob_i, 3];
		edge_i = edge_i + 1;
	end
	if i==n1,
		edge(edge_i,:) = [glob_i, 4];
		edge_i = edge_i + 1;
	end
	if j==1,
		edge(edge_i,:) = [glob_i, 1];
		edge_i = edge_i + 1;
	end
	if j==n2,
		edge(edge_i,:) = [glob_i, 2];
		edge_i = edge_i + 1;
	end
	if k==1,
		edge(edge_i,:) = [glob_i, 5];
		edge_i = edge_i + 1;
	end
	if k==n3,
		edge(edge_i,:) = [glob_i, 6];
		edge_i = edge_i + 1;
	end
	glob_i = glob_i + 1;
end
end
end

function [Ah] = stiffness_2D(dofs,p,tri)
% description:
%      generate the stiffness matrix Ah;
%
% arguments:
%   - dofs  Degrees of freedom.  
%   - p     nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   elements. Index to the three corners of element i given in row i.
% returns:
%		- Ah 	  Stiffness matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: March 2015

    Ah = sparse(dofs,dofs);
		%Ah = spalloc(dofs,dofs,6*dofs); %Estimate number of nnz elements, saves a lot of memory!
    % Nq = 4; %number of integration points in quadrature2D
    Ne = length(tri(:,1)); %number of elements
    
    for i = 1:Ne %iterates over all elements
        
        pis = tri(i,:);   %pi is node numbers in element i, f.ex. [7,8,1]
        p1 = p(pis(1),:); %coordinates to p1 in element i, f.ex. [3,1]
        p2 = p(pis(2),:); %coordinates to p2 in element i
        p3 = p(pis(3),:); %coordinates to p3 in element i
        
        X = [p1(1) p2(1) p3(1)];
        Y = [p1(2) p2(2) p3(2)];
        A_k = polyarea(X,Y); %Area of the element,(Jacobi determinant) 
        
        delPhi = zeros(2,3);
        delPhi(1,1) = p2(2)-p3(2);
        delPhi(2,1) = p3(1)-p2(1);
        delPhi(1,2) = p3(2)-p1(2);
        delPhi(2,2) = p1(1)-p3(1);
        delPhi(1,3) = p1(2)-p2(2);
        delPhi(2,3) = p2(1)-p1(1);
        
        delPhi = (1/(2*A_k))*delPhi;
        
        Ah(pis,pis) = Ah(pis,pis) + A_k*delPhi'*delPhi;

    end
end

function [Gh] = gradient(dofs,p,tri,b)
% description:
%      generate the gradient matrix Gh;
%
% arguments:
%   - dofs  Degrees of freedom.  
%   - p     Nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   Elements. Index to the three corners of element i given in row i.
%		- b 		function handle 2-dim Vector field 
% returns:
%		- Gh 	  Gradient matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: March 2015

    %Gh = sparse(dofs,dofs);
		Gh = spalloc(dofs,dofs,7*dofs); %Estimate number of nnz elements, saves a lot of memory!
    Nq = 4; %number of integration points in quadrature2D
    Ne = length(tri(:,1)); %number of elements
    
    for i = 1:Ne %iterates over all elements
        
        pis = tri(i,:);   %pi is node numbers in element i, f.ex. [7,8,1]
        p1 = p(pis(1),:); %coordinates to p1 in element i, f.ex. [3,1]
        p2 = p(pis(2),:); %coordinates to p2 in element i
        p3 = p(pis(3),:); %coordinates to p3 in element i
        
        X = [p1(1) p2(1) p3(1)];
        Y = [p1(2) p2(2) p3(2)];
        %A_k = 2*polyarea(X,Y); %Area of the element,(Jacobi determinant) 
        A_k = (p2(1)-p1(1))*(p3(2)-p1(2))-(p3(1)-p1(1))*(p2(2)-p1(2)); %Direct

        %Basis functions
        phi1 = @(x,y) (1/(A_k))*((p2(1)*p3(2) - p3(1)*p2(2)) + (p2(2)-p3(2))*x + (p3(1)-p2(1))*y);
        phi2 = @(x,y) (1/(A_k))*((p3(1)*p1(2) - p1(1)*p3(2)) + (p3(2)-p1(2))*x + (p1(1)-p3(1))*y);
        phi3 = @(x,y) (1/(A_k))*((p1(1)*p2(2) - p2(1)*p1(2)) + (p1(2)-p2(2))*x + (p2(1)-p1(1))*y); 
        PHI = @(x,y) [phi1(x,y) phi2(x,y) phi3(x,y)];
        
        % Derivatives of basis functions
        delPhi = zeros(2,3);
        delPhi(1,1) = p2(2)-p3(2);
        delPhi(2,1) = p3(1)-p2(1);
        delPhi(1,2) = p3(2)-p1(2);
        delPhi(2,2) = p1(1)-p3(1);
        delPhi(1,3) = p1(2)-p2(2);
        delPhi(2,3) = p2(1)-p1(1);
        delPhi = delPhi/A_k; 

        G = @(x,y) kron(b(x,y),PHI(x,y));
				I = quadrature2D(p1,p2,p3,Nq,G); 
        
        Gh(pis,pis) = Gh(pis,pis) + I'*delPhi;
    end
end



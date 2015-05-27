function [Dh] = gradient(dofs,p,tri,b,mu)
% description:
%      generate the gradient matrix Dh;
%
% arguments:
%   - dofs  Degrees of freedom.  
%   - p     nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   elements. Index to the three corners of element i given in row i.
%   - b     Vector field function handle 
%   - mu    Viscosity coefficient
% returns:
%		- Dh 	  "gradient" matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: March 2015

    %Dh = sparse(dofs,dofs);
		Dh = spalloc(dofs,dofs,14*dofs); %Estimate number of nnz elements, saves a lot of memory!
    Nq = 4; %number of integration points in quadrature2D
    Ne = length(tri(:,1)); %number of elements
    
    for i = 1:Ne %iterates over all elements
        
        Dh_loc = zeros(9);
        coords = [1 2 4 5 7 8];

        pis = tri(i,:);   %pi is node numbers in element i, f.ex. [7,8,1]
        p1 = p(pis(1),:); %coordinates to p1 in element i, f.ex. [3,1]
        p2 = p(pis(2),:); %coordinates to p2 in element i
        p3 = p(pis(3),:); %coordinates to p3 in element i
        
				% Local to Global mapping indices
				Map = reshape([3*(pis-1)+1;3*(pis-1)+2;3*pis],1,9); % Map is now a vector containing the corresponding global indices

        X = [p1(1) p2(1) p3(1)];
        Y = [p1(2) p2(2) p3(2)];
        %A_k = 2*polyarea(X,Y); % twice the area of the element,(Jacobi determinant) 
        A_k = abs((p1(1)-p3(1))*(p2(2)-p1(2)) - (p1(1)-p2(1))*(p3(2)-p1(2))); % Jacobian 

        %Basis functions
        phi1 = @(x,y) (1/(A_k))*((p2(1)*p3(2) - p3(1)*p2(2)) + (p2(2)-p3(2))*x + (p3(1)-p2(1))*y);
        phi2 = @(x,y) (1/(A_k))*((p3(1)*p1(2) - p1(1)*p3(2)) + (p3(2)-p1(2))*x + (p1(1)-p3(1))*y);
        phi3 = @(x,y) (1/(A_k))*((p1(1)*p2(2) - p2(1)*p1(2)) + (p1(2)-p2(2))*x + (p2(1)-p1(1))*y); 
        PHI = @(x,y) [phi1(x,y) phi2(x,y) phi3(x,y)];
        
        %Derivative of Basis functions
        delPhi = zeros(2,3);
        delPhi(1,1) = p2(2)-p3(2);
        delPhi(2,1) = p3(1)-p2(1);
        delPhi(1,2) = p3(2)-p1(2);
        delPhi(2,2) = p1(1)-p3(1);
        delPhi(1,3) = p1(2)-p2(2);
        delPhi(2,3) = p2(1)-p1(1);
        delPhi = delPhi/A_k;
        
        G = @(x,y) reshape(mu*delPhi-kron(b(x,y),PHI(x,y)),6,1);
        D = @(x,y) kron(G(x,y)',G(x,y));
				I = quadrature2D(p1,p2,p3,Nq,D); 
        Dh_loc(coords,coords) = I;
        Dh(Map,Map) = Dh(Map,Map) + Dh_loc;
    end
end



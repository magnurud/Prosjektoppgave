function [Dh] = stiffness_2D(dofs,p,tri,b,mu)
% description:
%      generate the gradient matrix Dh;
%
% arguments:
%   - dofs  Degrees of freedom.  
%   - p     nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   elements. Index to the three corners of element i given in row i.
%   - b     Vector field function handle (Assumed constant).
%   - mu    Viscosity coefficient
% returns:
%		- Dh 	  "gradient" matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: March 2015

    Dh = sparse(dofs,dofs);
		%Dh = spalloc(dofs,dofs,6*dofs); %Estimate number of nnz elements, saves a lot of memory!
    % Nq = 4; %number of integration points in quadrature2D
    Ne = length(tri(:,1)); %number of elements
    B = b(0,0);
    
    for i = 1:Ne %iterates over all elements
        
        pis = tri(i,:);   %pi is node numbers in element i, f.ex. [7,8,1]
        p1 = p(pis(1),:); %coordinates to p1 in element i, f.ex. [3,1]
        p2 = p(pis(2),:); %coordinates to p2 in element i
        p3 = p(pis(3),:); %coordinates to p3 in element i
        
				% Local to Global mapping indices
				Map = reshape([3*(pis-1)+1;3*(pis-1)+2;3*pis],1,9); % Map is now a vector containing the corresponding global indices

        X = [p1(1) p2(1) p3(1)];
        Y = [p1(2) p2(2) p3(2)];
        A_k = 2*polyarea(X,Y); % twice the area of the element,(Jacobi determinant) 
        
        delPhi = zeros(2,3);
        delPhi(1,1) = p2(2)-p3(2);
        delPhi(2,1) = p3(1)-p2(1);
        delPhi(1,2) = p3(2)-p1(2);
        delPhi(2,2) = p1(1)-p3(1);
        delPhi(1,3) = p1(2)-p2(2);
        delPhi(2,3) = p2(1)-p1(1);

        I1 = zeros(9,9); % This will be the stiffness element matrix for this element! 
        I2 = zeros(9,9); % This will be the stiffness element matrix for this element! 
				% to help out with numerating
				n1 = [1 4 7]; n2 = [2 5 8];	n3 = [3 6 9];

        %%%% D1 - first and easiest part of the matrix %%%%

				% Component 1,1   
				I1(n1,n1) = A_k/24*kron(ones(3,1)*B(1),ones(1,3)*B(1));
				% Component 1,2
				I1(n1,n2) = A_k/24*kron(ones(3,1)*B(1),ones(1,3)*B(2));
				% Component 2,1
				I1(n2,n1) = A_k/24*kron(ones(3,1)*B(2),ones(1,3)*B(1));
				% Component 2,2
				I1(n2,n2) = A_k/24*kron(ones(3,1)*B(2),ones(1,3)*B(2));

        % Correcting for the diagonal elements where i=j, doubling the
        % blockdiagonal elements
        I1 = I1 + diag(diag(I1,-1),-1);
        I1 = I1 + diag(diag(I1,0),0);
        I1 = I1 + diag(diag(I1,1),1);
        Dh(Map,Map) = Dh(Map,Map) + I1;  % Can do this directly as well.

        %%%% D2 - second and harder part of the matrix %%%%

				% Component 1,1   
				I2(n1,n1) = A_k/6*(kron(delPhi(1,:)',ones(1,3)*B(1))  +  kron(ones(3,1)*B(1),delPhi(1,:)));
				% Component 1,2
				I2(n1,n2) = A_k/6*(kron(delPhi(1,:)',ones(1,3)*B(2))  +  kron(ones(3,1)*B(1),delPhi(2,:)));
				% Component 2,1
				I2(n2,n1) = A_k/6*(kron(delPhi(2,:)',ones(1,3)*B(1))  +  kron(ones(3,1)*B(2),delPhi(1,:)));
				% Component 2,2
				I2(n2,n2) = A_k/6*(kron(delPhi(2,:)',ones(1,3)*B(2))  +  kron(ones(3,1)*B(2),delPhi(2,:)));

        Dh(Map,Map) = Dh(Map,Map) - mu*I2;  % Can do this directly as well.

        if(i == -11) 
          I2
          A_k*1000
          delPhi
        end


    end
end


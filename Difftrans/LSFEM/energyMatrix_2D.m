function [Ah] = energyMatrix_2D(dofs,p,tri)
% description:
%      generate the energy matrix Ah;
%
% arguments:
%   - dofs  Degrees of freedom.  
%   - p     nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   elements. Index to the three corners of element i given in row i.
% returns:
%		- Ah 	  Energy matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: March 2015

		Ah = spalloc(dofs,dofs,5*dofs); %Estimate number of nnz elements, saves a lot of memory!
    Ne = length(tri(:,1)); %number of elements

		% Quadrature variables begin 
		Nq = 2;
		w = [1,1];
		xq = [-1/sqrt(3);1/sqrt(3)]; 
		xq = 1/2*xq +1/2;
		% Quatrature variables end

    for i = 1:Ne %iterates over all elements
        
        pis = tri(i,:);   %pi is node numbers in element i, f.ex. [7,8,1]
        p1 = p(pis(1),:); %coordinates to p1 in element i, f.ex. [3,1]
        p2 = p(pis(2),:); %coordinates to p2 in element i
        p3 = p(pis(3),:); %coordinates to p3 in element i
        
        X = [p1(1) p2(1) p3(1)];
        Y = [p1(2) p2(2) p3(2)];
        A_k = polyarea(X,Y); %Area of the element,(Jacobi determinant) 
        % Making the D-Matrix
        D = zeros(3,3,Nq);
        D(2,1,:) = p2(2)-p3(2);
        D(3,1,:) = p3(1)-p2(1);
        D(2,2,:) = p3(2)-p1(2);
        D(3,2,:) = p1(1)-p3(1);
        D(2,3,:) = p1(2)-p2(2);
        D(3,3,:) = p2(1)-p1(1);
				%%%%%%%%%%%%%%%%%%%%%%%%%%% Rather try the other approach first %%%%%%%%%%%%%% 
				D(1,1,:) = 1-
        
        D = (1/(2*A_k))*D;
        
        Ah(pis,pis) = Ah(pis,pis) + A_k*D'*D;

    end
end


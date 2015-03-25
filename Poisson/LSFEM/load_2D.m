function [fh] = load_2D(dofs,p,tri,f)
% function [fh] = load_2D(dofs,p,tri,f),
% 
% description:
%     Evaluates the loading function f in all the nodes.
%
% arguments:
%		- f 		Loading function expression of two variables
%   - dofs  Degrees of freedom. 
%   - p     Nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   Elements. Index to the three corners of element i given in row i.
% returns:
%   - fh 		The loading function evaluated in all the nodes.
%
% author: Magnus Aa. Rud
% last edit: March 2015

    fh = zeros(dofs,1);
    Nq = 4;
    Ne = length(tri(:,1)); %Number of elements
    
    for i = 1:Ne
        pis = tri(i,:);
        p1 = p(pis(1),:);
        p2 = p(pis(2),:);
        p3 = p(pis(3),:);

        A_k = abs((p1(1)-p3(1))*(p2(2)-p1(2)) - (p1(1)-p2(1))*(p3(2)-p1(2))); % Jacobian

				Map = reshape([3*(pis-1)+1;3*(pis-1)+2;3*pis],1,9); % Map is now a vector containing the corresponding global indices
        
        delPhi = zeros(2,3);
        delPhi(1,1) = p2(2)-p3(2);
        delPhi(2,1) = p3(1)-p2(1);
        delPhi(1,2) = p3(2)-p1(2);
        delPhi(2,2) = p1(1)-p3(1);
        delPhi(1,3) = p1(2)-p2(2);
        delPhi(2,3) = p2(1)-p1(1);

        delPhi = (1/A_k)*delPhi;


				M = [delPhi;0 0 0];	
        
				I = quadrature2D(p1,p2,p3,Nq,f)
        % If f is const
        %I = A_k*f(0,0)/2;
				fh(Map) = fh(Map) + I*reshape(M,9,1);       
    end
end


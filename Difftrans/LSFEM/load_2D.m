function [fh] = load_2D(dofs,p,tri,f,b,mu)
% function [fh] = load_2D(dofs,p,tri,f,b,mu),
% 
% description:
%     Evaluates the loading function f in all the nodes.
%
% arguments:
%		- f 		Loading function expression of two variables
%   - dofs  Degrees of freedom. 
%   - p     Nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   Elements. Index to the three corners of element i given in row i.
%   - mu    Viscosity
%   - b     Vector field as a function handle
% returns:
%   - fh 		The loading function evaluated in all the nodes.
%
% author: Magnus Aa. Rud
% last edit: March 2015

    fh = zeros(dofs,1);
    Nq = 4;
    Ne = length(tri(:,1)); %Number of elements
    B = b; % Assuming constant vector field
    I = 0;
    
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

        delPhi = delPhi/A_k;

        % for a non-const f 
        %Basis functions
        phi1 = @(x,y) (1/(A_k))*((p2(1)*p3(2) - p3(1)*p2(2)) + (p2(2)-p3(2))*x + (p3(1)-p2(1))*y);
        phi2 = @(x,y) (1/(A_k))*((p3(1)*p1(2) - p1(1)*p3(2)) + (p3(2)-p1(2))*x + (p1(1)-p3(1))*y);
        phi3 = @(x,y) (1/(A_k))*((p1(1)*p2(2) - p2(1)*p1(2)) + (p1(2)-p2(2))*x + (p2(1)-p1(1))*y); 

        PHI = @(x,y) [phi1(x,y) phi2(x,y) phi3(x,y)];
        M = mu*[delPhi;0 0 0];
        F = @(x,y) (M-[B(1)*PHI(x,y) ; B(2)*PHI(x,y); 0 0 0])*f(x,y);
				I_2 = quadrature2D(p1,p2,p3,Nq,F); 
        fh(Map) = fh(Map) + reshape(I_2,9,1);       

        % If f is const
         %M = 1/2*mu*[delPhi;0 0 0];
         %M(1,:) = M(1,:) - ones(1,3)*B(1)/6;
         %M(2,:) = M(2,:) - ones(1,3)*B(2)/6;
         %I = A_k*f(0,0);
         %fh(Map) = fh(Map) + I*reshape(M,9,1);       
    end
end

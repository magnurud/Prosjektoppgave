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

        A_k = (1/2)*abs((p1(1)-p3(1))*(p2(2)-p1(2)) - (p1(1)-p2(1))*(p3(2)-p1(2)));
        
        %f=@(x,y) -8*pi*cos(2*pi*(x^2+y^2))+16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2)); %As argument in function
        
        for n = 1:3
            if n == 1
                phi = @(x,y) (1/(2*A_k))*((p2(1)*p3(2) - p3(1)*p2(2)) + (p2(2)-p3(2))*x + (p3(1)-p2(1))*y);
            elseif n == 2
                phi = @(x,y) (1/(2*A_k))*((p3(1)*p1(2) - p1(1)*p3(2)) + (p3(2)-p1(2))*x + (p1(1)-p3(1))*y);
            elseif n == 3
                phi = @(x,y) (1/(2*A_k))*((p1(1)*p2(2) - p2(1)*p1(2)) + (p1(2)-p2(2))*x + (p2(1)-p1(1))*y); 
            end   
            g = @(x,y) f(x,y)*phi(x,y);   
            I = quadrature2D(p1,p2,p3,Nq,g);
            fh(pis(n),1) = fh(pis(n),1) + I;       
        end         
    end
end


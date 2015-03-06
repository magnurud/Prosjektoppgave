%% This function calculates the elements in the stiffness matrix Ah
% [p tri edge] = getDisk(Nr), Nr = number of nodes
% tri = Nx3-matrix, N = number of elements, node numbers in columns
% p(i) = coordinates to node i
% edge = hold node numbers to edges

function [Ah] = stiffness_2D(Nr,p,tri)

    Ah = sparse(Nr,Nr);
    % Nq = 4; %number of integration points in quadrature2D
    N = length(tri(:,1)); %number of elements
    
    for i = 1:N %iterates over all elements
        
        pis = tri(i,:); %pi is node numbers in element i, f.eks. [7,8,1]
        p1 = p(pis(1),:); %coordinates to p1 in element i, f.eks. [3,1]
        p2 = p(pis(2),:); %coordinates to p2 in element i
        p3 = p(pis(3),:); %coordinates to p3 in element i
        
        X = [p1(1) p2(1) p3(1)];
        Y = [p1(2) p2(2) p3(2)];
        A_k = polyarea(X,Y); %Area of the element
        
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

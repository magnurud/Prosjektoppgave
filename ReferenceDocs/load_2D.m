function [fh] = load_2D(Nr,p,tri)
    fh = zeros(Nr,1);
    Nq = 4;
    N = length(tri(:,1));
    
    for i = 1:N
        pis = tri(i,:);
        p1 = p(pis(1),:);
        p2 = p(pis(2),:);
        p3 = p(pis(3),:);

        A_k = (1/2)*abs((p1(1)-p3(1))*(p2(2)-p1(2)) - (p1(1)-p2(1))*(p3(2)-p1(2)));
        
        f=@(x,y) -8*pi*cos(2*pi*(x^2+y^2))+16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2));
        
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
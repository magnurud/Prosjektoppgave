function I = quadrature2D(p1,p2,p3,Nq,g)
%g = @(x,y) log(x+y);
%p1 = [1;0]; p2 = [3;1]; p3 = [3;2];
j1 = p1-p3;
j2 = p2-p3;
Jacobi = [j1(1) j1(2); j2(1) j2(2)];
    if Nq == 1
        eta = [1/3 1/3 1/3];
        rho = 1;
    elseif Nq == 3
        eta = [0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5];
        rho = [1/3; 1/3; 1/3];
    elseif Nq == 4
        eta = [1/3 1/3 1/3; 3/5 1/5 1/5; 1/5 3/5 1/5; 1/5 1/5 3/5];
        rho = [-9/16; 25/48; 25/48; 25/48];
    end    
    I = 0;
    for i = 1:Nq
        x = eta(i,1)*p1 + eta(i,2)*p2 + eta(i,3)*p3;
        I = I + rho(i)*g(x(1),x(2));
    end
    I = 0.5*abs(det(Jacobi))*I; %1/2 since 2 dimensions and times the absolute value of Jacobi
    
end
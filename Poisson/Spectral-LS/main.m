% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2
% 		 using spectral-least squares methods;
%
% author: Magnus Aa. Rud
% last edit: April 2015

N = 20; %Number of points in each direction
NLS = 3*N; % Number of unknowns in each direction
dofs = N^2;
LSdofs = NLS^2;
f = @(x,y) 1; % Loading function
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
LDM = LagrangeDerivativeMatrix_GLL(N);
f = @(x,y) 5*pi^2*sin(pi*x)*sin(2*pi*y); % Loading function
u = @(x,y) sin(pi*x)*sin(2*pi*y); % Analytical solution
U = zeros(dofs,1);
for I = 1:dofs
  j = mod(I-1,N)+1;
  i = fix((I-1)/N)+1;
  U(I) = u(x(i),y(j));
end

% Assembling stiffness matrix

Ah = stiffness_2D(N,x,y,wX,wY,LDM);
fh = load_2D(N,x,y,wX,wY,f,LDM);

% Boundary conditions
for I = 1:dofs
	J = 3*I;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N)
    Ah(J,:) = 0;
    Ah(J,J) = 1;
    fh(J) = 1;
  elseif(j==1 || j==N)
    Ah(J,:) = 0;
    Ah(J,J) = 1;
    fh(J) = 1;
  end
end
    
uh = Ah\fh;


% Plotting
figure;
surf(x,y,reshape(uh(3:3:LSdofs)',N,N));
title('Numerical Solution');

figure;
surf(x,y,reshape(U,N,N));
title('Analytical Solution');

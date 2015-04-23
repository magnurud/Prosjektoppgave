% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2
% 		 using spectral-least squares methods;
%
% author: Magnus Aa. Rud
% last edit: April 2015

N = 30; %Number of points in each direction
h = 1/(N-1);
NLS = 3*N; % Number of unknowns in each direction
dofs = N^2;
LSdofs = 3*dofs;
f = @(x,y) 1; % Loading function
b = @(x,y) [1,0]; % vector field
mu = 10*1E-2;
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
W = diag(wX);
LDM = 2*LagrangeDerivativeMatrix_GLL(N);
B1 = zeros(N,N);
B2 = zeros(N,N);
for I = 1:dofs
  j = mod(I-1,N)+1;
  i = fix((I-1)/N)+1;
  bloc = b(x(i),y(j));
  B1(i,j) = bloc(1);
  B2(i,j) = bloc(2); 
end

Ah = stiffness_2D_fast(W,LDM,mu);
Gh = gradient_2D_fast(mu,B1,B2,W,LDM,dofs);
Ah = Ah+Gh;
fh = load_2D_fast(N,x,y,wX,wY,f,LDM);
% Boundary conditions
for I = 1:dofs
	J = I+2*dofs;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N || j==1 || j==N)
    Ah(J,:) = 0;
    Ah(:,J) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
  end
end

uh = Ah\fh;

% Plotting
figure;
surf(x,y,reshape(uh(2*dofs+1:end),N,N)');
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')


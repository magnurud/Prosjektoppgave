% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2 using spectral methods;
%
% author: Magnus Aa. Rud
% last edit: April 2015

N = 30; %Number of points in each direction
dofs = N^2;
f = @(x,y) 1; % Loading function
b = @(x,y) [1,0]; % vector field
mu = 1*1E-2;
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
LDM = 2*LagrangeDerivativeMatrix_GLL(N); % Need to multiply with 2/(b-a)
W = diag(wX);
f = @(x,y) 1 ;
B1 = zeros(N,N);
B2 = zeros(N,N);
for I = 1:dofs
  j = mod(I-1,N)+1;
  i = fix((I-1)/N)+1;
  bloc = b(x(i),y(j));
  B1(i,j) = bloc(1);
  B2(i,j) = bloc(2); 
end

% Assembling stiffness matrix

Ah = stiffness_2D_fast(wX,LDM);
Gh = gradient_2D(LDM,B1,B2,N,wX,wY);
Gh2 = gradient_2D_fast(LDM,B1,B2,W,dofs);
fh = load_2D(N,x,y,wX,wY,f);
Ah = mu*Ah+Gh2;

% Boundary conditions
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N || j==1 || j==N)
    Ah(I,:) = 0;
    Ah(I,I) = 1;
    fh(I) = 0;
  end
end
    
uh = Ah\fh;

% Plotting
figure;
surf(x,y,reshape(uh',N,N)');
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

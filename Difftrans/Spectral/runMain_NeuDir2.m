function [eh cn] = runMain_NeuDir2(N,mu,alpha)
% runMain.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2
% 		 using spectral-least squares methods;
%
% arguments:
%   - N     number of discretization points in each direction
%   - mu    the diffusion constant
%   - alpha The size of the vector field b 
%          
% returns:
%		- cn  Condition number for the matrix
%   - eh  Error compared to the amalytical solution
% author: Magnus Aa. Rud
% last edit: April 2015

dofs = N^2;
B = @(x,y) alpha*[2;1]; %Vector Field
f = @(x,y) 1;
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
LDM = 2*LagrangeDerivativeMatrix_GLL(N); % Need to multiply with 2/(b-a)
W = diag(wX);
B1 = zeros(N,N);
B2 = zeros(N,N);
U = zeros(dofs,1);
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  bloc = B(x(i),y(j));
  B1(i,j) = bloc(1);
  B2(i,j) = bloc(2); 
  U(I) = u(x(i),y(j));
end
Ah = stiffness_2D_fast(wX,LDM);
Gh = gradient_2D_fast(LDM,B1,B2,W,dofs);
fh = load_2D(N,x,y,wX,wY,f);
Ah = mu*Ah+Gh;

% Dirichlet
g_south = @(x,y) 1-x;
g_east  = @(x,y) y ;

% Neumann
h_west = @(x,y) 1;
h_north= @(x,y) 1;

% Boundary conditions
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if( j==1) 
		fh = fh - Ah(:,I)*g_south(x(i),y(j));
	end
	if( i==N)
		fh = fh - Ah(:,I)*g_east(x(i),y(j));
	end
	if( i==1)
		fh(I) = fh(I)+mu*wX(j)*h_west(x(i),y(j));
	end
	if( j==N)
		fh(I) = fh(I)+mu*wX(i)*h_north(x(i),y(j));
  end
end

for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(j==1) 
    Ah(I,:) = 0;
    Ah(:,I) = 0;
    Ah(I,I) = 1;
    fh(I) = g_south(x(i),y(j));
	end
	if(i==N)
    Ah(I,:) = 0;
    Ah(:,I) = 0;
    Ah(I,I) = 1;
    fh(I) = g_east(x(i),y(j));
  end
end
    
uh = Ah\fh;

% Plotting
figure;
surf(x,y,reshape(uh,N,N)');
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

cn = condest(Ah);
%eh = norm((uh-U),'inf')/norm(U,'inf');
eh = norm((uh-U))/norm(U);


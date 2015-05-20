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
%
% author: Magnus Aa. Rud
% last edit: April 2015

h = 1/(N-1);
NLS = 3*N; % Number of unknowns in each direction
dofs = N^2;
LSdofs = 3*dofs;
B = @(x,y) alpha*[1;0]; %Vector Field
f = @(x,y) 1; % loading function
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

Ah = stiffness_2D_fast(W,LDM,mu);
Gh = gradient_2D_fast(mu,B1,B2,W,LDM,dofs);
Ah = Ah+Gh;
fh = load_2D_fast(N,x,y,wX,wY,f,LDM,B1,B2,mu);

%% Dirchlet boundary conditions %%
g1 = @(x,y) 1-x; % South side boundary function
g2 = @(x,y) y; % East side boundary function

%% Neumann BC %% 
%The amount that leaves the domain
g3 = @(x,y) 1; % West side boundary function
g4 = @(x,y) 1; % North side boundary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Vectorized BC's %% 
Rg = zeros(LSdofs,1);
for I = 1:dofs
	J = I+2*dofs;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(j==1) 
    Rg(J) = g1(x(i),y(j)); % South side
  elseif(i == N) 
    Rg(J) = g2(x(i),y(j)); % East side
  elseif(i == 1) 
    Rg(I) = -g3(x(i),y(j)); % West side 
  elseif(j == N) 
    Rg(J-dofs) = g4(x(i),y(j)); % North side
  end
end

fh = fh - Ah*Rg;

% Boundary conditions
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==N || j==1)
		J = I+2*dofs;
    Ah(J,:) = 0;
    Ah(:,J) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
	elseif(i==1)
		J = I;
    Ah(J,:) = 0;
    Ah(:,J) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
	elseif(j==N)
		J = I+dofs;
    Ah(J,:) = 0;
    Ah(:,J) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;

  end
end

uh = Ah\fh;
uh = uh+Rg;

% Plotting
figure;
surf(x,y,reshape(uh(2*dofs+1:end),N,N)');
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

eh = norm((uh(2*dofs+1:end)-U),'inf')/norm(U,'inf');
cn = condest(Ah);



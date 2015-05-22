function [eh cn] = runMain(N,mu,alpha)
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
B = @(x,y) alpha*[1;1]; %Vector Field
% TEST 1 %
f = @(x,y) mu*exp(x)*(pi^2-1)*sin(pi*y)...
		+exp(x)*B(x,y)'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution

% TEST 2 %
%f = @(x,y) -mu*(6*x-pi^2*x^3)*sin(pi*y)...
    %+(B(x,y)')*[3*x^2*sin(pi*y) ; x^3*pi*cos(pi*y)]; % Loading function
%u = @(x,y) x^3*sin(pi*y); % Analytical solution
% TEST 3 %
%f = @(x,y) 1;
%u = @(x,y) 0; % Analytical solution
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
g1 = @(x,y) u(x,y); % South side boundary function
g2 = @(x,y) u(x,y); % East side boundary function
g3 = @(x,y) u(x,y); % West side boundary function
g4 = @(x,y) u(x,y); % North side boundary function
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
    Rg(J) = g3(x(i),y(j)); % West side
  elseif(j == N) 
    Rg(J) = g4(x(i),y(j)); % North side
  end
end
fh = fh - Ah*Rg;

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
uh = uh+Rg;

%% Plotting
%figure;
%subplot(1,2,1)
%surf(x,y,reshape(uh(2*dofs+1:end),N,N)');
%title('Numerical Solution');
%xlabel('x')
%ylabel('y')
%zlabel('z')

%%% Plotting the analytical solution
%subplot(1,2,2) % second subplot
%surf(x,y,reshape(U,N,N)');
%title('Analytical Solution');
%xlabel('x')
%ylabel('y')
%zlabel('z')

eh = norm((uh(2*dofs+1:end)-U),'inf')/norm(U,'inf');
cn = condest(Ah);
%Peclet number
Peclet = alpha*h/(2*mu);

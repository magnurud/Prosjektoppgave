function [eh cn] = runMain(N)
% runMain.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2
% 		 using spectral-least squares methods;
%
% author: Magnus Aa. Rud
% last edit: April 2015

h = 1/(N-1);
NLS = 3*N; % Number of unknowns in each direction
dofs = N^2;
LSdofs = 3*dofs;
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
W = diag(wX);
LDM = 2*LagrangeDerivativeMatrix_GLL(N);
f = @(x,y) exp(x)*(pi^2-1)*sin(pi*y);
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
U = zeros(dofs,1);
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  U(I) = u(x(i),y(j));
end

% Assembling stiffness matrix
% The variables are structured block-wise
Ah = stiffness_2D_fast(W,LDM);
fh = load_2D_precond(N,x,y,wX,wY,f,LDM);

%% Dirchlet boundary conditions %%
g1 = @(x,y) u(x,y); % South side boundary function
g2 = @(x,y) u(x,y); % East side boundary function
g3 = @(x,y) u(x,y); % West side boundary function
g4 = @(x,y) u(x,y); % North side boundary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Vectorized BC's %% 
Rg = zeros(LSdofs,1);
for I = 1:dofs
  j = fix((I-1)/N)+1;
  i = mod(I-1,N)+1;
  if(j==1) 
    Rg(I+2*dofs) = g1(x(i),y(j)); % South side
  elseif(i == N) 
    Rg(I+2*dofs) = g2(x(i),y(j)); % East side
  elseif(i == 1) 
    Rg(I+2*dofs) = g3(x(i),y(j)); % West side
  elseif(j == N) 
    Rg(I+2*dofs) = g4(x(i),y(j)); % North side
  end
end

fh = fh - Ah*Rg;

% Boundary conditions
for I = 1:dofs
	J = I+2*dofs;
  j = mod(I-1,N)+1;
  i = fix((I-1)/N)+1;
  if(i==1 || i==N || j==1 || j==N)
    Ah(J,:) = 0;
    Ah(:,J) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
  end
end

uh = Ah\fh;
uh = uh + Rg;

eh = norm(uh(2*dofs+1:end)-U,'inf')/norm(U,'inf');
norm(uh(2*dofs+1:end)-U,'inf')
norm(U,'inf')
cn = cond(Ah);

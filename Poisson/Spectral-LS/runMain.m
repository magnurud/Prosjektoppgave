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
f = @(x,y) 1; % Loading function
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
W = diag(wX);
LDM = 2*LagrangeDerivativeMatrix_GLL(N);
f = @(x,y) 5*pi^2*sin(pi*x)*sin(2*pi*y); % Loading function
u = @(x,y) sin(pi*x)*sin(2*pi*y); % Analytical solution
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


% Boundary conditions
for I = 1:dofs
	J = I+2*dofs;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N || j==1 || j==N)
    Ah(J,:) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
  end
end

uh = Ah\fh;

eh = norm(uh(2*dofs+1:end)-U,'inf')/norm(U,'inf');
cn = cond(Ah);

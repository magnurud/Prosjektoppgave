function cn = runMain(N,mu,alpha)
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

delta = 0.1;
h = 1/(N-1);
NLS = 3*N; % Number of unknowns in each direction
dofs = N^2;
LSdofs = 3*dofs;
f = @(x,y) 1; % Loading function
b = @(x,y) alpha*[1,0]; % vector field
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

% Assembling LS part %
A_LS = stiffness_LS(W,LDM,mu);
G_LS = gradient_LS(mu,B1,B2,W,LDM,dofs);
A_LS = A_LS+G_LS;
f_LS = load_LS(N,x,y,wX,wY,f,LDM);

% Assembling Spectral part 
A_Sp = stiffness_Spec(W,LDM);
G_Sp = gradient_Spec(LDM,B1,B2,W,dofs);
f_Sp = load_Spec(N,x,y,wX,wY,f);
A_Sp = mu*Ah+Gh2;

%%%%% COMBINING %%%%%
A_LS = delta*A_LS;
f_LS = delta*f_LS;

A_LS(2*dofs+1:end,2*dofs+1:end) = A_LS(2*dofs+1:end,2*dofs+1:end) + A_Sp;
f_LS(2*dofs+1:end) = f_LS(2*dofs+1:end) + f_Sp;
%%%%% DONE COMBINING %%%%%

Ah = sparse(A_LS);
fh = sparse(f_LS);
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
cn = condest(Ah);

% Plotting
figure;
surf(x,y,reshape(uh(2*dofs+1:end),N,N)');
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

%Peclet number
Peclet = max(max(sqrt(B1.^2+B2.^2)))*h/(2*mu)

function [eh cn] = runMain(N,mu,alpha,sigma,delta)
% runMain.m
%
% description:
%      Solving the diff-trans-reaction problem on the square (0,1)^2
% 		 using spectral-GLS method;
%
% arguments:
%   - N     number of discretization points in each direction
%   - mu    the diffusion constant
%   - alpha The size of the vector field b 
%   - sigma the weight of the reaction term
%   - delta the weight of the spectral part
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
f = @(x,y) 1;
u = @(x,y) 0;% exp(x)*sin(pi*y); % Analytical solution
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

% Assembling LS part %
A_LS = stiffness_LS(W,LDM,mu);
G_LS = gradient_LS(mu,B1,B2,W,LDM,dofs);
R_LS = reaction_LS(mu,B1,B2,sigma,W,LDM,dofs);
A_LS = A_LS+G_LS+R_LS;
f_LS = load_LS(N,x,y,wX,wY,f,LDM,B1,B2,mu,sigma);

% Assembling Spectral part 
A_Sp = stiffness_Spec(W,LDM);
G_Sp = gradient_Spec(LDM,B1,B2,W,dofs);
R_Sp = reaction_Spec(W);
f_Sp = load_Spec(N,x,y,wX,wY,f);
A_Sp = mu*A_Sp+G_Sp+sigma*R_Sp;

%%%%% COMBINING %%%%%
A_LS = delta*A_LS;
f_LS = delta*f_LS;

A_LS(2*dofs+1:end,2*dofs+1:end) = A_LS(2*dofs+1:end,2*dofs+1:end) + A_Sp;
f_LS(2*dofs+1:end) = f_LS(2*dofs+1:end) + f_Sp;
%%%%% DONE COMBINING %%%%%

Ah = sparse(A_LS);
fh = sparse(f_LS);

% Homogenous Boundary conditions
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

eh = norm((uh(2*dofs+1:end)-U),'inf')/norm(U,'inf');
cn = condest(Ah);

%Peclet number
Peclet = max(max(sqrt(B1.^2+B2.^2)))*h/(2*mu)


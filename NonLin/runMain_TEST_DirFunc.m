function [eh cn] = runMain_TEST_DirFunc(N,mu,alpha)
% runMain.m
%
% description:
% 		 Testing the functions! 
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
tic
maxit = 15;
eVec = zeros(maxit+1,1);
Convrate = zeros(maxit,1);
eVec(1)=1;
h = 1/(N-1);
dofs = N^2;
LSdofs = 3*dofs;
sigma = 0;
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
w2 = @(x,y) pi*exp(x)*cos(pi*y);
dB = @(x,y) alpha*[x;y]; %Vector Field, linear part
f = @(x,y) exp(x)*(mu*pi^2-mu)*sin(pi*y)...
    +exp(x)*(dB(x,y)*u(x,y))'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square

LDM = 2*LagrangeDerivativeMatrix_GLL(N); % Need to multiply with 2/(b-a)
W = diag(wX);
dB1 = zeros(dofs); dB2 = zeros(dofs);
BB1 = zeros(dofs); BB2 = zeros(dofs);
U = zeros(dofs,1); %Analytical solution
U_anal = zeros(LSdofs,1); % Total analytical solution
uh = 1*ones(LSdofs,1); % Initial guess
F  = zeros(dofs,1); % Loading function 
for I = 1:dofs
	i = mod(I-1,N)+1;
	j = fix((I-1)/N)+1;
	% uh(2*dofs+I) = 0.3*u(x(i),y(j));
	U(I) = u(x(i),y(j));
	U_anal(I) = -u(x(i),y(j));
	U_anal(dofs+I) = -w2(x(i),y(j));
	U_anal(2*dofs+I) = u(x(i),y(j));
	% The constant dB-matrices
	bloc = dB(x(i),y(j));
	dB1(I,I) = bloc(1);
	dB2(I,I) = bloc(2); 
	BB1(I,I) = bloc(1)*U(I);
	BB2(I,I) = bloc(2)*U(I);
	F(I) = f(x(i),y(j));
end
%uh = 0.5*U_anal;

% Assembling the linear matrices
A_LS = stiffness_LS(W,LDM,mu);
A_L = sparse(A_LS);
F_L = zeros(LSdofs,1);
%% Dirchlet boundary conditions %%
g1 = @(x,y) u(x,y); % South side boundary function
g2 = @(x,y) u(x,y); % East side boundary function
g3 = @(x,y) u(x,y); % West side boundary function
g4 = @(x,y) u(x,y); % North side boundary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for I = 1:dofs
	J = I+2*dofs;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(j==1) 
    A_L(J,J) = A_L(J,J) + W(i,i); % South side
		F_L(J) = F_L(J) + W(i,i)*g1(x(i),y(j)); 
  elseif(i == N) 
    A_L(J,J) = A_L(J,J) + W(j,j); % East side
		F_L(J) = F_L(J) + W(j,j)*g2(x(i),y(j)); 
  elseif(i == 1) 
    A_L(J,J) = A_L(J,J) + W(j,j); % West side
		F_L(J) = F_L(J) + W(j,j)*g3(x(i),y(j)); 
  elseif(j == N) 
    A_L(J,J) = A_L(J,J) + W(i,i); % North side
		F_L(J) = F_L(J) + W(i,i)*g4(x(i),y(j)); 
  end
end

	U1 = U_anal(2*dofs+1:end);
	% fictive Updating the vectorfield matrix
	B1 = diag(dB1*(U1)); 
	B2 = diag(dB2*(U1));

	% Getting the u-dependent matrices
	G_LS = gradient_LS(mu,B1,B2,W,LDM,dofs);
	F_NL = load_LS2(W,LDM,B1,B2,F,mu,sigma);
	A_NL = sparse(G_LS); % The non-linear part of A, is already 

	% Calculating the final u-dependent term corresponding to the loading function.
	fh = F_L+F_NL;%+A_NL*[zeros(2*dofs,1);Rg];

  uh = (A_NL+A_L)\fh;
  eh = norm((uh(2*dofs+1:end)-U),'inf')/norm(U,'inf');
  %eh = norm(uh+[zeros(2*dofs,1) ; Rg ] - U_anal,'inf');%/norm(U_anal,'inf');

% Plotting
if(1)
figure;
subplot(1,2,1)
%surf(x,y,reshape(uh,N,N)');
surf(x,y,reshape(uh(2*dofs+1:end),N,N)');
%surf(x,y,reshape(uh(1:dofs),N,N)');
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

%% Plotting the analytical solution
subplot(1,2,2) % second subplot
surf(x,y,reshape(U,N,N)');
title('Analytical Solution');
xlabel('x')
ylabel('y')
zlabel('z')
end
cn = condest(A_L+A_NL);

%Peclet number
%Peclet = max(max(sqrt(B1.^2+B2.^2)))*h/(2*mu)
toc






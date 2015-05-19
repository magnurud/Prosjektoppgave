function [J eh cn] = runMain_LS(N,mu,alpha)
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
%   - sigma the weight of the reaction term
%   - delta the weight of the spectral part
% returns:
%		- cn  Condition number for the matrix
%
% author: Magnus Aa. Rud
% last edit: April 2015
tic
maxit = 18;
eVec = zeros(maxit+1,1);
Convrate = zeros(maxit,1);
eVec(1)=1;
h = 1/(N-1);
dofs = N^2;
LSdofs = 3*dofs;
sigma = 0;
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
w2 = @(x,y) pi*exp(x)*cos(pi*y);
dB = @(x,y) alpha*[x^2;1-2*y]; %Vector Field, linear part
f = @(x,y) exp(x)*(mu*pi^2-mu)*sin(pi*y)...
    +exp(x)*(dB(x,y)*u(x,y))'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square

LDM = 2*LagrangeDerivativeMatrix_GLL(N); % Need to multiply with 2/(b-a)
W = diag(wX);
dB1 = zeros(dofs);
dB2 = zeros(dofs);
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
F(I) = f(x(i),y(j));
end
uh = 0.5*U_anal;
%% Dirchlet boundary conditions %%
g1 = @(x,y) u(x,y); % South side boundary function
g2 = @(x,y) u(x,y); % East side boundary function
g3 = @(x,y) u(x,y); % West side boundary function
g4 = @(x,y) u(x,y); % North side boundary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Vectorized BC's %% 
Rg = zeros(dofs,1);
for I = 1:dofs
	J = I;
	i = mod(I-1,N)+1;
	j = fix((I-1)/N)+1;
	if(j==1) 
		Rg(J) = g1(x(i),y(j)); % South side
		uh(J+2*dofs) = 0;
	elseif(i == N) 
		Rg(J) = g2(x(i),y(j)); % East side
		uh(J+2*dofs) = 0;
	elseif(i == 1) 
		Rg(J) = g3(x(i),y(j)); % West side
		uh(J+2*dofs) = 0;
	elseif(j == N) 
		Rg(J) = g4(x(i),y(j)); % North side
		uh(J+2*dofs) = 0;
	end
end

% Assembling the linear matrices
A_LS = stiffness_LS(W,LDM,mu);
A_L = sparse(A_LS);
F_L = A_L*[zeros(2*dofs,1);Rg];

% Homogenous Boundary conditions
for I = 1:dofs
	J = 2*dofs+I;
	i = mod(I-1,N)+1;
	j = fix((I-1)/N)+1;
	if(i==1 || i==N || j==1 || j==N)
		A_L(J,:) = 0;
		A_L(:,J) = 0;
		A_L(J,J) = 1;
	end
end

% Code to repeat
for it = 1:maxit
	U1 = uh(2*dofs+1:end);
	W1 = uh(1:dofs);
	W2 = uh(dofs+1:2*dofs);

	% Updating the vectorfield matrix
	B1 = diag(dB1*(U1+Rg)); 
	B2 = diag(dB2*(U1+Rg));

	% Getting the u-dependent matrices
	%G_Sp = gradient_Spec(LDM,B1,B2,W,dofs);
	%F_NL = load_LS(N,x,y,wX,wY,f,LDM,B1,B2,mu,sigma);
	G_LS = gradient_LS(mu,B1,B2,W,LDM,dofs);
	F_NL = load_LS2(W,LDM,B1,B2,F,mu,sigma);
	A_NL = sparse(G_LS); % The non-linear part of A

	% Getting the components of the jacobi matrices
	%J_Sp = jacobi_Spec_lin(W,LDM,dB1,dB2,B1,B2,U1,Rg);
	J_LS = jacobi_LS_lin(W,LDM,dB1,dB2,B1,B2,W1,W2,mu);
	J_f  = jacobi_LS_load_lin(W,dB1,dB2,F);
  J_LS = sparse(J_LS-J_f); 

	% Calculating the final u-dependent term corresponding to the loading function.
	fh = F_L-F_NL;%+A_NL*[zeros(2*dofs,1);Rg];

  % Homogenous Boundary conditions
	% The one non-zero term in each row is taken care of by A_L
  for I = 1:dofs
    J = 2*dofs+I;
    i = mod(I-1,N)+1;
    j = fix((I-1)/N)+1;
    if(i==1 || i==N || j==1 || j==N)
      A_NL(J,:) = 0;
      A_NL(:,J) = 0;
      J_LS(J,:) = 0;
      J_LS(:,J) = 0;
      %J_LS(J,J) = 1;
      fh(J) = 0;
    end
  end
	% Step 1
	r = (A_L+A_NL)*uh+fh;
	%Step 2
	J = sparse(A_L+A_NL+J_LS); % The J-matrix used in the newtons iteration
	e = J\r;
	%Step 3
	uh = uh-e;
	%U1 = uh;
	% Boundary adjustment
  uh_BC = uh(2*dofs+1:end)+Rg;
  eh = norm((uh_BC-U),'inf')/norm(U,'inf');
  %eh = norm(uh+[zeros(2*dofs,1) ; Rg ] - U_anal,'inf');%/norm(U_anal,'inf');
  eVec(it+1) = eh;
  Convrate(it) = eh/((eVec(it)^2));
  
end
uh = uh_BC;
format long
Convrate

figure(1);
plot(1:maxit,log(eVec(2:end)))

% Plotting
if(0)
figure;
subplot(1,2,1)
surf(x,y,reshape(uh,N,N)');
%surf(x,y,reshape(uh(dofs+1:2*dofs),N,N)');
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
eh = norm((uh-U),'inf')/norm(U,'inf');
cn = condest(A_L+A_NL);

%Peclet number
%Peclet = max(max(sqrt(B1.^2+B2.^2)))*h/(2*mu)
toc



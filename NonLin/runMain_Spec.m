function [eh cn] = runMain_Spec(N,mu,alpha,sigma,delta)
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
h = 1/(N-1);
dofs = N^2;
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
dB = @(x,y) alpha*[2*x;-y]; %Vector Field, linear part
f = @(x,y) exp(x)*(mu*pi^2-mu)*sin(pi*y)...
    +exp(x)*(dB(x,y)*u(x,y))'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square

LDM = 2*LagrangeDerivativeMatrix_GLL(N); % Need to multiply with 2/(b-a)
W = diag(wX);
dB1 = zeros(dofs);
dB2 = zeros(dofs);
U = zeros(dofs,1); %Analytical solution
uh = ones(dofs,1); % Initial guess
U1 = uh;
F  = zeros(dofs,1); % Loading function 
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  uh(I) = 0.3*u(x(i),y(j));
  U(I) = u(x(i),y(j));
% The constant dB-matrices
bloc = dB(x(i),y(j));
dB1(I,I) = bloc(1);
dB2(I,I) = bloc(2); 
F(I) = f(x(i),y(j));
end

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
	elseif(i == N) 
		Rg(J) = g2(x(i),y(j)); % East side
	elseif(i == 1) 
		Rg(J) = g3(x(i),y(j)); % West side
	elseif(j == N) 
		Rg(J) = g4(x(i),y(j)); % North side
	end
end

% Assembling the linear matrices
A_Sp = mu*stiffness_Spec(W,LDM);
A_L = sparse(A_Sp);

%f_Sp = load_Spec(N,x,y,wX,wY,f);
f_Sp = load_Spec2(W,F);
% The gradient part due to non-homogenous BC's
	B1R = diag(dB1*Rg); 
	B2R = diag(dB2*Rg);
	G_R = gradient_Spec(LDM,B1R,B2R,W,dofs);

%F_L = A_L*Rg + G_R*Rg -f_Sp;
F_L = A_L*Rg -f_Sp;

% Homogenous Boundary conditions
for I = 1:dofs
	J = I;
	i = mod(I-1,N)+1;
	j = fix((I-1)/N)+1;
	if(i==1 || i==N || j==1 || j==N)
		A_L(J,:) = 0;
		A_L(:,J) = 0;
		A_L(J,J) = 1;
		F_L(J) = 0;
	end
end

% Code to repeat
for it = 1:15 
	% Updating the vectorfield matrix
	B1 = diag(dB1*(U1+Rg)); 
	B2 = diag(dB2*(U1+Rg));

	G_Sp = gradient_Spec(LDM,B1,B2,W,dofs);
	J_Sp = jacobi_Spec_lin(W,LDM,dB1,dB2,B1,B2,U1,Rg);
	A_NL = sparse(G_Sp); % The non-linear part of A, is already 
  	J_LS = sparse(J_Sp);
	fh 	 = F_L+A_NL*Rg;

  % Homogenous Boundary conditions
	% The one non-zero term in each row is taken care of by A_L
  for I = 1:dofs
    J = I;
    i = mod(I-1,N)+1;
    j = fix((I-1)/N)+1;
    if(i==1 || i==N || j==1 || j==N)
      A_NL(J,:) = 0;
      A_NL(:,J) = 0;
      J_LS(J,:) = 0;
      J_LS(:,J) = 0;
      J_LS(J,J) = 1;
      fh(J) = 0;
    end
  end
	% Step 1
	r = (A_L+A_NL)*uh+fh;
	%Step 2
	B = sparse(A_L+J_LS); % The B-matrix used in the newtons iteration
	e = B\r;
	%Step 3
	uh = uh-e;
	U1 = uh;
	% Boundary adjustment
  uh_BC = uh+Rg;
  eh = norm((uh_BC-U),'inf')/norm(U,'inf');
end
uh = uh_BC;


% Plotting
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

eh = norm((uh-U),'inf')/norm(U,'inf');
cn = condest(A_L+A_NL);

%Peclet number
%Peclet = max(max(sqrt(B1.^2+B2.^2)))*h/(2*mu)
toc


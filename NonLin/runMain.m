function [eh cn] = runMain(N,mu,alpha,sigma,delta)
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
NLS = 3*N; % Number of unknowns in each direction
dofs = N^2;
LSdofs = 3*dofs;
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
dB = @(x,y) alpha*[2*x;-y]; %Vector Field, linear part
f = @(x,y) exp(x)*(mu*pi^2-mu)*sin(pi*y)...
    +exp(x)*(dB(x,y)*u(x,y))'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
LDM = 2*LagrangeDerivativeMatrix_GLL(N); % Need to multiply with 2/(b-a)
W = diag(wX);
dB1 = zeros(N,N);
dB2 = zeros(N,N);
U = zeros(dofs,1); %Analytical solution
uh = ones(dofs,1); % Initial guess
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  uh(I) = 0.3*u(x(i),y(j));
  U(I) = u(x(i),y(j));
	% The constant dB-matrices
	bloc = dB(x(i),y(j));
	dB1(i,j) = bloc(1);
	dB2(i,j) = bloc(2); 
end

% Assembling the linear matrices
A_LS = stiffness_LS(W,LDM,mu);
A_Sp = mu*stiffness_Spec(W,LDM);
A_LS = delta*A_LS;
A_LS(2*dofs+1:end,2*dofs+1:end) = A_LS(2*dofs+1:end,2*dofs+1:end) + A_Sp;
A_L = sparse(A_LS);

f_Sp = load_Spec(N,x,y,wX,wY,f);
% Homogenous Boundary conditions
for I = 1:dofs
	J = I+2*dofs;
	i = mod(I-1,N)+1;
	j = fix((I-1)/N)+1;
	if(i==1 || i==N || j==1 || j==N)
		A_L(J,:) = 0;
		A_L(:,J) = 0;
		A_L(J,J) = 1;
		f_Sp(J) = 0;
	end
end
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

% Code to repeat
for it = 1:15 
	% Updating the vectorfield matrix
	B1 = diag(dB1*U1); 
	B2 = diag(dB2*U1);
%% THE NON-LINEAR PART %%
	G_LS = gradient_LS(mu,B1,B2,W,LDM,dofs);
	G_Sp = gradient_Spec(LDM,B1,B2,W,dofs);
  f_LS = load_LS(N,x,y,wX,wY,f,LDM,B1,B2,mu,sigma);
%% THE JACOBIAN PART %%
  J_LS  = jacobi_LS_lin(W,LDM,dB1,dB2,B1,B2,U1,W1,mu);
	J_Sp  = jacobi_Spec_lin(W,LDM,dB1,dB2,B1,B2,U1);
  Jf_LS = jacobi_LS_load_lin(N,x,y,wX,wY,f,LDM,dB1,dB2);
	
  %%%%% COMBINING NON-LIN PART %%%%%
  G_LS = delta*G_LS;
  f_LS = delta*f_LS;

  G_LS(2*dofs+1:end,2*dofs+1:end) = G_LS(2*dofs+1:end,2*dofs+1:end) + G_Sp;
  f_LS(2*dofs+1:end) = f_LS(2*dofs+1:end) + f_Sp;
	A_NL = sparse(G_LS); % The non-linear part of A
	fh = sparse(f_LS); 

  %%%%% COMBINING JACOBI PART %%%%%
  J_LS = delta*(Jf_LS-J_LS);
  J_LS(2*dofs+1:end,2*dofs+1:end) = J_LS(2*dofs+1:end,2*dofs+1:end) - J_Sp;
  %%%%% DONE COMBINING %%%%%

	%  fh = fh - Ah*Rg;

  % Homogenous Boundary conditions
  for I = 1:dofs
    J = I+2*dofs;
    i = mod(I-1,N)+1;
    j = fix((I-1)/N)+1;
    if(i==1 || i==N || j==1 || j==N)
      A_NL(J,:) = 0;
      A_NL(:,J) = 0;
      A_NL(J,J) = 1;
      fh(J) = 0;
    end
  end
	% Step 1
	r = fh-A_L*uh-A_NL*uh;
	%Step 2
	B = sparse(A_L-J_LS); % The B-matrix used in the newtons iteration
	e = B\r;
	%Step 3
	uh = uh+e;
	U1 = uh(2*dofs+1:end);
	W1 = uh(1:dofs);
	% Boundary adjustment
  %uh = uh+Rg;
  eh = norm((uh(2*dofs+1:end)-U),'inf')/norm(U,'inf')
end


% Plotting
figure;
subplot(1,2,1)
surf(x,y,reshape(uh(2*dofs+1:end),N,N)');
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

eh = norm((uh(2*dofs+1:end)-U),'inf')/norm(U,'inf');
cn = condest(Ah);

%Peclet number
Peclet = max(max(sqrt(B1.^2+B2.^2)))*h/(2*mu)
toc

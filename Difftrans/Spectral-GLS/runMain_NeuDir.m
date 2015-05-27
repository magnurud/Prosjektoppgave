function [eh cn] = runMain_NeuDir(N,mu,alpha,delta)

	% runMain.m
	%
	% description:
	%      Solving the Diffusion transport problem on the square (0,1)^2
	% 		 using spectral-least squares methods;
	%
	% arguments:
	%   - N     number of discretization points in each direction
	%   - mu    the diffusion constant
	%   - alpha The size of the vector field b 
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
	f = @(x,y) mu*exp(x)*(pi^2-1)*sin(pi*y)...
	+exp(x)*B(x,y)'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
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

	% Assembling LS part %
	A_LS = stiffness_LS(W,LDM,mu);
	G_LS = gradient_LS(mu,B1,B2,W,LDM,dofs);
	A_LS = A_LS+G_LS;
	f_LS = load_LS(N,x,y,wX,wY,f,LDM,B1,B2,mu);

	% Assembling Spectral part 
	A_Sp = stiffness_Spec(W,LDM);
	G_Sp = gradient_Spec(LDM,B1,B2,W,dofs);
	f_Sp = load_Spec(N,x,y,wX,wY,f);
	A_Sp = mu*A_Sp+G_Sp;

	%%%%% COMBINING %%%%%
	A_LS = delta*A_LS;
	f_LS = delta*f_LS;

	A_LS(2*dofs+1:end,2*dofs+1:end) = A_LS(2*dofs+1:end,2*dofs+1:end) + A_Sp;
	f_LS(2*dofs+1:end) = f_LS(2*dofs+1:end) + f_Sp;
	%%%%% DONE COMBINING %%%%%

	Ah = sparse(A_LS);
	fh = sparse(f_LS);
	%%Dirichlet boundary conditions %%
	g1 = @(x,y) u(x,y); % South side boundary function
	g2 = @(x,y) u(x,y); % East side boundary function
	g3 = @(x,y) u(x,y); % West side boundary function
	%% Neumann boundary conditions %%
	g4 = @(x,y) pi*exp(x)*cos(pi*y) ; % North side boundary function
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	%% Vectorized BC's %% 
	Rg = zeros(LSdofs,1);
	for I = 1:dofs
		i = mod(I-1,N)+1;
		j = fix((I-1)/N)+1;
		J = 2*dofs+I;
		if(j==1) 
			Rg(J) = g1(x(i),y(j)); % South side
		elseif(i == N) 
			Rg(J) = g2(x(i),y(j)); % East side
		elseif(i == 1) 
			Rg(J) = g3(x(i),y(j)); % West side
		% Neumann part
		elseif(j == N) 
			Rg(J-dofs) = -g4(x(i),y(j)); % North side
			fh(J) = fh(J)+1/mu*wX(i)*g4(x(i),y(j));
		end
  end
		%%%%%%%%%%%%% 
	fh = fh - Ah*Rg;

	% Homogenous Neumann Boundary conditions
	for I = 1:dofs
		i = mod(I-1,N)+1;
		j = fix((I-1)/N)+1;
		if(i==1 || i==N || j==1) % Dirichlet on the solution vector
			J = 2*dofs+I;
			Ah(J,:) = 0;
			Ah(:,J) = 0;
			Ah(J,J) = 1;
			fh(J) = 0;
		end
		if(j==N) % Dirichlet on the y-comp of w 
			J = I+dofs;
			Ah(J,:) = 0;
			Ah(:,J) = 0;
			Ah(J,J) = 1;
			fh(J) = 0;
			%%%%%%%%%% 
		end
	end

	uh = Ah\fh;
	uh = uh+Rg;

	% Plotting
	figure;
	subplot(1,2,1)
	surf(x,y,reshape(uh(2*dofs+1:end),N,N)');
	%surf(x,y,reshape(uh(dofs+1:2*dofs),N,N)');
	%surf(x,y,reshape(uh(1:dofs),N,N)');
	%surf(x,y,reshape(Rg(1:dofs),N,N)');
	%surf(x,y,reshape(Rg(dofs+1:2*dofs),N,N)');
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


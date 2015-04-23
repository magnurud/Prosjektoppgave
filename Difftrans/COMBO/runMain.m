function cn = runMain(N,mu,alpha)
% runMain.m
%
% description:
%      Solving the diffusion transport problem on the square (0,1)^2;
%
% arguments:
%   - N     number of grid points in each direction
%   - mu    diffusion coeffisient
%   - alpha vector field constant
% returns:
%		- cn    Condition number of the stiffness matrix
%
% author: Magnus Aa. Rud
% last edit: March 2015

dofs = 3*N^2; % Number of degrees of freedom.
NN = N^2; %Number of nodes
[p,tri,e] = getSquare(N); %nodes, edges and elements.
f = @(x,y) 1; % Loading function
b = @(x,y) alpha*[1 ; 1]; % Vector field creating the transport

% Assemble Matrices and loading function for LSFEM 
	Ah = stiffness_2D(dofs,p,tri,mu); % Now needs to include viscosity
	fh_LS= load_2D(dofs,p,tri,f,b,mu);
	Dh = gradient_2D(dofs,p,tri,b,mu);
	K_LS = Ah+Dh; % Total matrix
	%K = Ah; % Total matrix
% 

% Assemble Matrices and loading function for FEM
	Ah = stiffness_2D_FEM(dofs/3,p,tri);
	fh_FEM = load_2D_FEM(dofs/3,p,tri,f);
	Gh = gradient_2D_FEM(dofs/3,p,tri,b(0,0));
	K_FEM = mu*Ah+Gh; % Total matrix
% 


% Imposing Dirchlet homogenous boundary conditions
for i = e
  K_FEM(i,:) = 0;
  K_FEM(:,i) = 0;
  K_FEM(i,i) = 1;
  fh_FEM(i) = 0;
end

% Imposing Homogenous boundary conditions 
for j = e
	i = 3*j;	
  K_LS(i,:) = 0;
  K_LS(:,i) = 0;
  K_LS(i,i) = 1;
  fh_LS(i) = 0;
end

%
K_LS(3:3:dofs,3:3:dofs) = K_LS(3:3:dofs,3:3:dofs) + K_FEM;
fh_LS(3:3:dofs) = fh_LS(3:3:dofs) + fh_FEM;
% Solving 
K = K_LS;
fh = fh_LS;
uh = K\fh;
	
%% Plotting the numerical solution
figure(1);
trisurf(tri,p(:,1),p(:,2),uh(3:3:dofs));
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

uh_max = uh(3:3:dofs);
cn = condest(K);
%norm(uh(3:3:dofs)-U)/norm(U);


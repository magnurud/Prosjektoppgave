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
Nodes = N^2;
NN = N^2; %Number of nodes
[p,tri,e] = getSquare(N); %nodes, edges and elements.
%f = @(x,y) 1; % Loading function
f = @(x,y) 5*pi^2*sin(pi*x)*sin(2*pi*y); % Loading function
b = @(x,y) alpha*[1 ; 1]; % Vector field creating the transport

u = @(x,y) sin(pi*x)*sin(2*pi*y); % Analytical solution
U = zeros(Nodes,1);
for I = 1:Nodes
    U(I) = u(p(I,1),p(I,2));
end

% Assemble Matrices and loading function for LSFEM 
	Ah = stiffness_2D(dofs,p,tri,mu); % Now needs to include viscosity
	fh= load_2D(dofs,p,tri,f,b,mu);
	Dh = gradient_2D(dofs,p,tri,b,mu);
	K = Ah+Dh; % Total matrix
	%K = Ah; % Total matrix
% 

% Imposing Homogenous boundary conditions 
for j = e
	i = 3*j;	
  K(i,:) = 0;
  K(:,i) = 0;
  K(i,i) = 1;
  fh(i) = 0;
end

%
% Solving 
uh = K\fh;
	
%% Plotting the numerical solution
figure(1);
trisurf(tri,p(:,1),p(:,2),uh(3:3:dofs));
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

%% Plotting the analytical solution
%figure(2);
%trisurf(tri,p(:,1),p(:,2),U);
%title('analytical Solution');
%xlabel('x')
%ylabel('y')
%zlabel('z')

uh_max = uh(3:3:dofs);
cn = condest(K);
%norm(uh(3:3:dofs)-U)/norm(U);


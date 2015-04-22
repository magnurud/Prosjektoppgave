function uh_max = runMain(N,mu)
% function uh = runMain(N,mu)
% runMain.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2 using LSFEM;
%
% author: Magnus Aa. Rud
% last edit: March 2015

dofs = 3*N^2; % Number of degrees of freedom.
NN = N^2; %Number of nodes
[p,tri,e] = getSquare(N); %nodes, edges and elements.
f = @(x,y) 1; % Loading function
b = @(x,y) [1 ; 0]; % Vector field creating the transport
%b = @(x,y) 0*[1 ; 1]; % Vector field creating the transport

% Assemble Matrices and loading function
	Ah = stiffness_2D(dofs,p,tri,mu); % Now needs to include viscosity
	fh = load_2D(dofs,p,tri,f,b,mu);
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

uh_max = uh(3:3:dofs);
%norm(uh(3:3:dofs)-U)/norm(U);


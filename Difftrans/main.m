% main.m
%
% description:
%      Solving the diffusion transport problem on the square (0,1)^2;
%
% author: Magnus Aa. Rud
% last edit: March 2015

N = 40; % Number of Nodes in each direction.
dofs = N^2; % Number of degrees of freedom.
[p,tri,e] = getSquare(N); %nodes, edges and elements.
f = @(x,y) 1; % Loading function
b = [-1 ; 1]; % Vector field creating the transport
mu = 1E-2; %Viscosity
%f=@(x,y) -8*pi*cos(2*pi*(x^2+y^2))+16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2)); %Loading function

% Assemble Matrices and loading function
	Ah = stiffness_2D(dofs,p,tri);
	fh = load_2D(dofs,p,tri,f);
	Gh = gradient_2D(dofs,p,tri,b);
	K = mu*Ah+Gh; % Total matrix
% 


% Imposing Dirchlet homogenous boundary conditions
for i = e
  K(i,:) = 0;
  K(:,i) = 0;
  K(i,i) = 1;
  fh(i) = 0;
end

%
%

% Solving 
uh = K\fh;
	
%% Plotting the numerical solution
figure(1);
trisurf(tri,p(:,1),p(:,2),uh);
title('Numerical Solution');


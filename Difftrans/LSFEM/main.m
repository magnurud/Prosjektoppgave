% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2 using LSFEM;
%
% author: Magnus Aa. Rud
% last edit: March 2015

N = 30; % Number of Nodes in each direction.
dofs = 3*N^2; % Number of degrees of freedom.
[p,tri,e] = getSquare(N); %nodes, edges and elements.
f = @(x,y) 1; % Loading function
%f=@(x,y) -8*pi*cos(2*pi*(x^2+y^2))+16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2)); %Loading function
b = @(x,y) [1 ; 1]; % Vector field creating the transport
mu = 0.0001; %Viscosity

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


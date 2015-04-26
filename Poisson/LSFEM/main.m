% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2 using LSFEM;
%
% author: Magnus Aa. Rud
% last edit: March 2015

N = 28; % Number of Nodes in each direction.
dofs = 3*N^2; % Number of degrees of freedom.
[p,tri,e] = getSquare(N); %nodes, edges and elements.
NN = length(e); %Number of boundary elements
f = @(x,y) 1; % Loading function
%f=@(x,y) -8*pi*cos(2*pi*(x^2+y^2))+16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2)); %Loading function

% Assemble Matrices and loading function
	Ah = stiffness_2D(dofs,p,tri);
	fh = load_2D(dofs,p,tri,f);
% 

% Imposing Homogenous boundary conditions the FEM Way
for j = e
	i = 3*j;	
  Ah(i,:) = 0;
  Ah(:,i) = 0;
  Ah(i,i) = 1;
  fh(i) = 0;
end

% Solving 
uh = Ah\fh;
	
%% Plotting the numerical solution
figure(1);
trisurf(tri,p(:,1),p(:,2),uh(3:3:dofs));
title('Numerical Solution');


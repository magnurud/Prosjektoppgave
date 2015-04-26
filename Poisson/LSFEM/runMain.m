function [eh cn] = runMain(N)
% function uh_max = runMain(N)
% description:
%      Solving the Poisson problem on the square (0,1)^2 with N discretization
%      points in each direction;
%
% author: Magnus Aa. Rud
% last edit: March 2015

dofs = 3*N^2; % Number of degrees of freedom.
Nodes = N^2;
[p,tri,e] = getSquare(N); %nodes, edges and elements.
f = @(x,y) 5*pi^2*sin(pi*x)*sin(2*pi*y); % Loading function
u = @(x,y) sin(pi*x)*sin(2*pi*y); % Analytical solution
U = zeros(Nodes,1);
for I = 1:Nodes
    U(I) = u(p(I,1),p(I,2));
end
% Assemble Matrices and loading function
	Ah = stiffness_2D(dofs,p,tri);
	fh = load_2D(dofs,p,tri,f);
% Imposing Homogenous boundary conditions
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
% figure(1);
% trisurf(tri,p(:,1),p(:,2),uh(3:3:dofs));
% title('Numerical Solution');
% figure(2);
% trisurf(tri,p(:,1),p(:,2),U);
% title('Analytical Solution');
eh = norm((uh(3:3:dofs)-U),'inf')/norm(U,'inf');
cn = condest(Ah);

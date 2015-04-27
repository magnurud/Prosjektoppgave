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

NN = 4*(N-1);
f = @(x,y) exp(x)*(pi^2-1)*sin(pi*y);
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
U = zeros(Nodes,1);
for I = 1:Nodes
    U(I) = u(p(I,1),p(I,2));
end
% Assemble Matrices and loading function
	Ah = stiffness_2D(dofs,p,tri);
	fh = load_2D(dofs,p,tri,f);

%% Dirchlet boundary conditions %%
g1 = @(x,y) u(x,y); % South side boundary function
g2 = @(x,y) u(x,y); % East side boundary function
g3 = @(x,y) u(x,y); % West side boundary function
g4 = @(x,y) u(x,y); % North side boundary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Vectorized BC's %% 
Rg = zeros(dofs,1);
for i = 1:NN
  E = fix(i/N); % Which edge
  j = e(i);
  if(E==0) 
    Rg(3*j) = g1(p(j,1),p(j,2));
  elseif(E==1) 
    Rg(3*j) = g2(p(j,1),p(j,2));
  elseif(E==2) 
    Rg(3*j) = g3(p(j,1),p(j,2));
  elseif(E==3) 
    Rg(3*j) = g4(p(j,1),p(j,2));
  end
end

fh = fh - Ah*Rg;

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
uh = uh+Rg;
%% Plotting the numerical solution
 figure(1);
 trisurf(tri,p(:,1),p(:,2),uh(3:3:dofs));
 title('Numerical Solution');
 figure(2);
 trisurf(tri,p(:,1),p(:,2),U);
 title('Analytical Solution');
eh = norm((uh(3:3:dofs)-U),'inf')/norm(U,'inf');
cn = condest(Ah);

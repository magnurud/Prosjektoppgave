% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2;
%
% author: Magnus Aa. Rud
% last edit: March 2015

N = 20; % Number of Nodes in each direction.
NN = 4*(N-1);
dofs = N^2; % Number of degrees of freedom.
[p,tri,e] = getSquare(N); %nodes, edges and elements.
%f = @(x,y) 5*pi^2*sin(pi*x)*sin(2*pi*y); % Loading function
f = @(x,y) exp(x)*(pi^2-1)*sin(pi*y);
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
U = zeros(dofs,1);
for I = 1:dofs
    U(I) = u(p(I,1),p(I,2));
end

% Assemble Matrices and loading function
	Ah = stiffness_2D(dofs,p,tri);
	fh = load_2D(dofs,p,tri,f);
% 

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
    Rg(j) = g1(p(j,1),p(j,2));
  elseif(E==1) 
    Rg(j) = g2(p(j,1),p(j,2));
  elseif(E==2) 
    Rg(j) = g3(p(j,1),p(j,2));
  elseif(E==3) 
    Rg(j) = g4(p(j,1),p(j,2));
  end
end

fh = fh - Ah*Rg;


% Imposing Homogenous boundary conditions
for i = e
  Ah(i,:) = 0;
  Ah(:,i) = 0;
  Ah(i,i) = 1;
  fh(i) = 0;
end

%
%

% Solving 
uh = Ah\fh;
uh = uh+Rg;
	
%% Plotting the numerical solution
figure(1);
trisurf(tri,p(:,1),p(:,2),uh);
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')


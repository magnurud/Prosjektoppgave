function [eh cn] = runMain(N,mu,alpha)
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

dofs = N^2; % Number of degrees of freedom.
NN = 4*(N-1);
[p,tri,e] = getSquare(N); %nodes, edges and elements.
b = alpha*[1 ; 1]; % Vector field creating the transport
B = @(x,y) alpha[1;1];
f = @(x,y) mu*exp(x)*(pi^2-1)*sin(pi*y)...
    +exp(x)*(b(1)*sin(pi*y)+pi*b(2)*cos(pi*y)); % Loading function
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
U = zeros(dofs,1);
for I = 1:dofs
    U(I) = u(p(I,1),p(I,2));
end

% Assemble Matrices and loading function
	Ah = stiffness_2D(dofs,p,tri);
	fh = load_2D(dofs,p,tri,f);
	Gh = gradient_2D(dofs,p,tri,b);
	Gh2 = gradient(dofs,p,tri,B);
	K = mu*Ah+Gh; % Total matrix
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

% Imposing Dirchlet homogenous boundary conditions to the non-homogenized system
for i = e
  K(i,:) = 0;
  K(:,i) = 0;
  K(i,i) = 1;
  fh(i) = 0;
end
%
% Solving 
uh = K\fh;
uh = uh+Rg;
	
%% Plotting the numerical solution
figure(1);
trisurf(tri,p(:,1),p(:,2),uh);
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

%% Plotting the analytical solution
figure(2);
trisurf(tri,p(:,1),p(:,2),U);
title('Analytical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

eh = norm((uh-U),'inf')/norm(U,'inf');
cn = condest(K);

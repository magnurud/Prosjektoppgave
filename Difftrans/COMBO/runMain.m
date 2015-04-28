function [eh cn] = runMain(N,mu,alpha,delta)
% runMain.m
%
% description:
%      Solving the diffusion transport problem on the square (0,1)^2;
%
% arguments:
%   - N     number of grid points in each direction
%   - mu    diffusion coeffisient
%   - alpha vector field constant
%   - delta the weight of the smoothing LS part
% returns:
%		- cn    Condition number of the stiffness matrix
%
% author: Magnus Aa. Rud
% last edit: March 2015

h = 1/(N-1);
Nodes = N^2;
dofs = 3*N^2; % Number of degrees of freedom.
NN = 4*(N-1); %Number of nodes
[p,tri,e] = getSquare(N); %nodes, edges and elements.
B = @(x,y) alpha*[x ; y]; % Vector field creating the transport
f = @(x,y) mu*exp(x)*(pi^2-1)*sin(pi*y)...
    +exp(x)*B(x,y)'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
%f = @(x,y) 1;
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
U = zeros(Nodes,1);
for I = 1:Nodes
    U(I) = u(p(I,1),p(I,2));
end

% Assemble Matrices and loading function for LSFEM 
	Ah = stiffness(dofs,p,tri,mu); % Now needs to include viscosity
	Dh = gradient(dofs,p,tri,B,mu);
	fh_LS= load_2D(dofs,p,tri,f,B,mu);
	K_LS = Ah+Dh; % Total matrix
% 

% Assemble Matrices and loading function for FEM
	Ah = stiffness_2D_FEM(dofs/3,p,tri);
	Gh = gradient_2D_FEM(dofs/3,p,tri,B);
	fh_FEM = load_2D_FEM(dofs/3,p,tri,f);
	K_FEM = mu*Ah+Gh; % Total matrix
% 

% Scaling the amount of smoothing
K_LS = delta*K_LS; fh_LS = delta*fh_LS;

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

  fh_LS = fh_LS - K_LS*Rg;
  fh_FEM= fh_FEM - K_FEM*Rg(3:3:dofs);
% Is now solving the corresponding homogenous problem

K_LS(3:3:dofs,3:3:dofs) = K_LS(3:3:dofs,3:3:dofs) + K_FEM;
fh_LS(3:3:dofs) = fh_LS(3:3:dofs) + fh_FEM;
K = K_LS; fh = fh_LS;

% Imposing Dirchlet homogenous boundary conditions
for i = e
	j = 3*i;	
  K(j,:) = 0;
  K(:,j) = 0;
  K(j,j) = 1;
  fh(j) = 0;
end

uh = K\fh;
uh = uh+Rg;

% Plotting
figure;
subplot(1,2,1)
trisurf(tri,p(:,1),p(:,2),uh(3:3:dofs));
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

%% Plotting the analytical solution
subplot(1,2,2) % second subplot
trisurf(tri,p(:,1),p(:,2),U);
title('Analytical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

	
%% Plotting the numerical solution
figure(1);

cn = condest(K);
eh = norm(uh(3:3:dofs)-U)/norm(U);
%eh = norm((uh(2*dofs+1:end)-U),'inf')/norm(U,'inf');
Peclet = norm(alpha*h/(2*mu))

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

dofs = 3*N^2; % Number of degrees of freedom.
Nodes = N^2;
NN = 4*(N-1); %Number of nodes
[p,tri,e] = getSquare(N); %nodes, edges and elements.
b = alpha*[1; 1]; % Vector field creating the transport
B = @(x,y) alpha*[x; y]; % Vector field creating the transport
f = @(x,y) mu*exp(x)*(pi^2-1)*sin(pi*y)...
    +exp(x)*B(x,y)'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
U = zeros(Nodes,1);
for I = 1:Nodes
    U(I) = u(p(I,1),p(I,2));
end

% Assemble Matrices and loading function for LSFEM 
	%Ah = stiffness_2D(dofs,p,tri,mu); % Now needs to include viscosity
	%Dh = gradient_2D(dofs,p,tri,b,mu);
    Ah2 = stiffness(dofs,p,tri,mu); % Now needs to include viscosity
    Dh2 = gradient(dofs,p,tri,B,mu);
    K = Ah2+Dh2; % Total matrix
    %K = Ah; % Total matrix

    fh= load_2D(dofs,p,tri,f,B,mu);
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
      Rg(3*j) = g1(p(j,1),p(j,2));
    elseif(E==1) 
      Rg(3*j) = g2(p(j,1),p(j,2));
    elseif(E==2) 
      Rg(3*j) = g3(p(j,1),p(j,2));
    elseif(E==3) 
      Rg(3*j) = g4(p(j,1),p(j,2));
    end
  end

  fh = fh - K*Rg;

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
  uh = uh+Rg;
    
  %%% Plotting the numerical solution
  %figure(1);
  %trisurf(tri,p(:,1),p(:,2),uh(3:3:dofs));
%title('Numerical Solution');
%xlabel('x')
%ylabel('y')
%zlabel('z')

% Plotting the analytical solution
%figure(2);
%trisurf(tri,p(:,1),p(:,2),U);
%title('analytical Solution');
%xlabel('x')
%ylabel('y')
%zlabel('z')

uh_max = uh(3:3:dofs);
eh = norm((uh(3:3:dofs)-U),'inf')/norm(U,'inf');
cn = condest(K);
%norm(uh(3:3:dofs)-U)/norm(U);


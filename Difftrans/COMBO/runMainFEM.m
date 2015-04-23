function cn = runMain(N,mu,alpha)
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
[p,tri,e] = getSquare(N); %nodes, edges and elements.
f = @(x,y) 1; % Loading function
b = alpha*[1 ; 0]; % Vector field creating the transport

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
xlabel('x')
ylabel('y')
zlabel('z')

%% Plotting the analytical solution
%figure(2);
%trisurf(tri,p(:,1),p(:,2),U);
%title('Analytical Solution');
cn = condest(K);

uh_max = uh;
xlabel('x')
ylabel('y')
zlabel('z')

function [cn] = runPlot(N,mu,alpha,delta,name)
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
f = @(x,y) 1;
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

cn = condest(K);
Peclet = norm(alpha*h/(2*mu))

figure;
trisurf(tri,p(:,1),p(:,2),uh(3:3:dofs));
xlabel('x')
ylabel('y')
zlabel('z')
view(40.5,30);
fig = gcf;

Filename = ['../../Latex/Figures/' name]
set(fig, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig, Filename , 'pdf') %Save figure


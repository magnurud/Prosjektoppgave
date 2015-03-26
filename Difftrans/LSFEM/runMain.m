function uh_max = runMain(N,mu)
% function uh = runMain(N,mu)
% runMain.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2 using LSFEM;
%
% author: Magnus Aa. Rud
% last edit: March 2015

dofs = 3*N^2; % Number of degrees of freedom.
NN = N^2; %Number of nodes
[p,tri,e] = getSquare(N); %nodes, edges and elements.
b = @(x,y) [1 ; 1]; % Vector field creating the transport
b(1,2)
b(2)
%f = @(x,y) 1; % Loading function
f = @(x,y) mu*5*pi^2*sin(pi*x)*sin(2*pi*y)...
          +b(1)*pi*cos(pi*x)*sin(2*pi*y)...
          +b(2)*2*pi*sin(pi*x)*cos(2*pi*y);
          f(0,0)
u = @(x,y) sin(pi*x)*sin(2*pi*y);
U = zeros(NN,1);
for I = 1:NN
    U(I) = u(p(I,1),p(I,2));
end

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

%% Plotting the analytical solution
figure(2);
trisurf(tri,p(:,1),p(:,2),U);
title('Analytical Solution');

uh_max = norm(uh(3:3:dofs)-U)/norm(U);


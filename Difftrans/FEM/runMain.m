function uh_max = runMain(N,mu)
% runMain.m
%
% description:
%      Solving the diffusion transport problem on the square (0,1)^2;
%
% author: Magnus Aa. Rud
% last edit: March 2015

dofs = N^2; % Number of degrees of freedom.
[p,tri,e] = getSquare(N); %nodes, edges and elements.
% PROBLEM // THIS GIVES THE RIGHT SOLN // %%%%%% 
%b = [1 ; -1]; % Vector field creating the transport
b = [1 ; 1];
%f = @(x,y) 1; % Loading function
f = @(x,y) mu*5*pi^2*sin(pi*x)*sin(2*pi*y)...
          +b(1)*pi*cos(pi*x)*sin(2*pi*y)...
          +b(2)*2*pi*sin(pi*x)*cos(2*pi*y);
u = @(x,y) sin(pi*x)*sin(2*pi*y);
U = zeros(dofs,1);
for I = 1:dofs
    U(I) = u(p(I,1),p(I,2));
end
%U = reshape(reshape(U,N,N)',dofs,1);
%f=@(x,y) -8*pi*cos(2*pi*(x^2+y^2))+16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2)); %Loading function

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
%figure(1);
%trisurf(tri,p(:,1),p(:,2),uh);
%title('Numerical Solution');

%% Plotting the analytical solution
%figure(2);
%trisurf(tri,p(:,1),p(:,2),U);
%title('Analytical Solution');

uh_max = norm(uh-U)/norm(U);

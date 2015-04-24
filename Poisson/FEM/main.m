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
f = @(x,y) 5*pi^2*sin(pi*x)*sin(2*pi*y); % Loading function
%f = @(x,y) 3;
%f = @(x,y) cos(2*pi*x)*y;
%f=@(x,y) -8*pi*cos(2*pi*(x^2+y^2))+16*pi^2*(x^2+y^2)*sin(2*pi*(x^2+y^2)); %Loading function

% Assemble Matrices and loading function
	Ah = stiffness_2D(dofs,p,tri);
	fh = load_2D(dofs,p,tri,f);
% 

%% Dirchlet boundary conditions %%
g1 = @(x,y) 0; % South side boundary function
g2 = @(x,y) y; % East side boundary function
g3 = @(x,y) y; % West side boundary function
g4 = @(x,y) 1; % North side boundary function
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


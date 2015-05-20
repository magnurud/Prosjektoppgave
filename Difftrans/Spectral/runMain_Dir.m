function [eh cn] = runMain_Dir(N,mu,alpha)
% runMain.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2
% 		 using spectral-least squares methods;
%
% arguments:
%   - N     number of discretization points in each direction
%   - mu    the diffusion constant
%   - alpha The size of the vector field b 
%          
% returns:
%		- cn  Condition number for the matrix
%   - eh  Error compared to the amalytical solution
% author: Magnus Aa. Rud
% last edit: April 2015

dofs = N^2;
B = @(x,y) alpha*[2;1]; %Vector Field
f = @(x,y) mu*exp(x)*(pi^2-1)*sin(pi*y)...
    +exp(x)*B(x,y)'*[sin(pi*y) ; pi*cos(pi*y)]; % Loading function
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
LDM = 2*LagrangeDerivativeMatrix_GLL(N); % Need to multiply with 2/(b-a)
W = diag(wX);
B1 = zeros(N,N);
B2 = zeros(N,N);
U = zeros(dofs,1);
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  bloc = B(x(i),y(j));
  B1(i,j) = bloc(1);
  B2(i,j) = bloc(2); 
  U(I) = u(x(i),y(j));
end
Ah = stiffness_2D_fast(wX,LDM);
Gh = gradient_2D_fast(LDM,B1,B2,W,dofs);
fh = load_2D(N,x,y,wX,wY,f);
Ah = mu*Ah+Gh;

% Boundary conditions
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N || j==1 || j==N)
		fh = fh - Ah(:,I)*u(x(i),y(j));
    Ah(I,:) = 0;
    Ah(:,I) = 0;
    Ah(I,I) = 1;
    fh(I) = u(x(i),y(j));
  end
end
    
uh = Ah\fh;

% Plotting
figure;
subplot(1,2,1) % first subplot
surf(x,y,reshape(uh,N,N)');
title('Numerical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

%% Plotting the analytical solution
subplot(1,2,2) % second subplot
surf(x,y,reshape(U,N,N)');
title('Analytical Solution');
xlabel('x')
ylabel('y')
zlabel('z')

cn = condest(Ah);
%eh = norm((uh-U),'inf')/norm(U,'inf');
eh = norm((uh-U))/norm(U);


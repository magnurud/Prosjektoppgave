function [eh cn] = runMain(N)
% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2 using spectral methods;
%
% author: Magnus Aa. Rud
% last edit: April 2015

dofs = N^2;
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
LDM = 2*LagrangeDerivativeMatrix_GLL(N); % Need to multiply with 2/(b-a)
f = @(x,y) exp(x)*(pi^2-1)*sin(pi*y);
u = @(x,y) exp(x)*sin(pi*y); % Analytical solution
U = zeros(dofs,1);
for I = 1:dofs
  j = mod(I-1,N)+1;
  i = fix((I-1)/N)+1;
  U(I) = u(x(i),y(j));
end

% Assembling stiffness matrix

Ah = stiffness_2D(N,x,y,wX,wY,LDM);
fh = load_2D(N,x,y,wX,wY,f);

%% Dirchlet boundary conditions %%
g1 = @(x,y) u(x,y); % South side boundary function
g2 = @(x,y) u(x,y); % East side boundary function
g3 = @(x,y) u(x,y); % West side boundary function
g4 = @(x,y) u(x,y); % North side boundary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Vectorized BC's %% 
Rg = zeros(dofs,1);
for I = 1:dofs
  i = fix((I-1)/N)+1;
  j = mod(I-1,N)+1;
  if(j==1) 
    Rg(I) = g1(x(i),y(j)); % South side
  elseif(i == N) 
    Rg(I) = g2(x(i),y(j)); % East side
  elseif(i == 1) 
    Rg(I) = g3(x(i),y(j)); % West side
  elseif(j == N) 
    Rg(I) = g4(x(i),y(j)); % North side
  end
end

fh = fh - Ah*Rg;

% Boundary conditions
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N || j==1 || j==N)
    Ah(I,:) = 0;
    Ah(:,I) = 0;
    Ah(I,I) = 1;
    fh(I) = 0;
  end
end
    
uh = Ah\fh;
uh = uh + Rg;
eh = norm(uh-U,'inf')/norm(U,'inf');
cn = cond(Ah);


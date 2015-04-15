% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2
% 		 using spectral-least squares methods;
%
% author: Magnus Aa. Rud
% last edit: April 2015

N = 10; %Number of points in each direction
NLS = 3*N; % Number of unknowns in each direction
dofs = N^2;
LSdofs = 3*dofs;
f = @(x,y) 1; % Loading function
[x,wX] = GLL_(N,0,1); % getting the GLL-points for the unit square
[y,wY] = GLL_(N,0,1); % getting the GLL-points for the unit square
W = diag(wX);
LDM = 2*LagrangeDerivativeMatrix_GLL(N);
f = @(x,y) 5*pi^2*sin(pi*x)*sin(2*pi*y); % Loading function
u = @(x,y) sin(pi*x)*sin(2*pi*y); % Analytical solution
U = zeros(dofs,1);
for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  U(I) = u(x(i),y(j));
end

% Assembling stiffness matrix

% The variables are structured node-wise
Ah = stiffness_2D(N,x,y,wX,wY,LDM);
fh = load_2D(N,x,y,wX,wY,f,LDM);

% The variables are structured block-wise
Ah2 = stiffness_2D_fast(W,LDM);
fh2 = load_2D_precond(N,x,y,wX,wY,f,LDM);

%%% The components of the block-diagonal matrix calculated in stif_fast;
%%% Only included here for test-purposes
A11 = kron(W,LDM'*W*LDM)+kron(W,W);

A12 = kron((W*LDM),(W*LDM)');

A13 = kron(W,W*LDM);

A22 = kron(LDM'*W*LDM,W)+kron(W,W);

A23 = kron(W*LDM,W);

A33 = kron(LDM'*W*LDM,W)+kron(W,LDM'*W*LDM);

%%% A preconditioner ... %%%
Ploc = kron(inv(W),eye(N));
P = blkdiag(Ploc,Ploc,Ploc);
Ah2 = P*Ah2;
fh2 = P*fh2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Boundary conditions
for I = 1:dofs
	J = 3*I;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N)
    Ah(J,:) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
  elseif(j==1 || j==N)
    Ah(J,:) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
  end
end
    
% Boundary conditions
for I = 1:dofs
	J = I+2*dofs;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N)
    Ah2(J,:) = 0;
    Ah2(J,J) = 1;
    fh2(J) = 0;
  elseif(j==1 || j==N)
    Ah2(J,:) = 0;
    Ah2(J,J) = 1;
    fh2(J) = 0;
  end
end
uh = Ah\fh;
uh2 = Ah2\fh2;




% Plotting
%figure;
%surf(x,y,reshape(uh(3:3:LSdofs),N,N));
%title('Numerical Solution');

%figure;
%surf(x,y,reshape(U,N,N));
%title('Analytical Solution');

error1 = norm(uh(3:3:LSdofs)-U,'inf')/norm(U,'inf')
error2 = norm(uh2(2*dofs+1:LSdofs)-U,'inf')/norm(U,'inf')
conditionNumber1 = cond(Ah)
conditionNumber2 = cond(Ah2)

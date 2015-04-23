% main.m
%
% description:
%      Solving the Poisson problem on the square (0,1)^2
% 		 using spectral-least squares methods;
%
% author: Magnus Aa. Rud
% last edit: April 2015

N = 15; %Number of points in each direction
h = 1/(N-1);
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

alpha = 10;
W_sq = sqrt(W);
Wt_sq = diag(1./sqrt(diag(W)));
I = eye(N);
%phi = inv(LDM'*LDM'+I);
phi = inv(LDM'*LDM'+alpha*I);
% Test of phi !  HERE IS THE BIG PROBLEM
'test of phi'
max(max(abs(phi*(LDM'*LDM'+alpha*I)-eye(N))))

% FIRST - CHOL DECOMP OF A 
L11 = kron(W_sq,complex(LDM'*W_sq,W_sq));
L22 = kron(W_sq,W_sq);
L33 = kron(I,I);
L21 = kron(LDM'*W_sq,W_sq);
L31 = kron(W_sq,complex(0,LDM'*W_sq));
L32 = kron(LDM'*W_sq,W_sq);
ZERO = zeros(dofs);
ONE = eye(dofs);
L = [L11 ZERO ZERO; L21 L22 ZERO ; L31 L32 ONE]; % Cholesky decomp of A! 
A = real(L*L');
B = imag(L*L');

% NEXT - THE INVERSE %
Lt11 = kron(Wt_sq,Wt_sq*phi*complex(LDM',-I));
Lt22 = kron(Wt_sq,Wt_sq);
Lt33 = kron(I,I);
Lt21 = -kron(Wt_sq*LDM',Wt_sq*phi*complex(LDM',-I));
Lt31 = kron(LDM'*LDM',phi*complex(LDM',-I)) - kron(I,LDM'*phi*complex(I,LDM'));
Lt32 = -kron(LDM',I);
ZERO = zeros(dofs);

Lt = [Lt11 ZERO ZERO; Lt21 Lt22 ZERO ; Lt31 Lt32 ONE]; % inverse of L!  
Lt2 = [Lt11 ZERO ZERO; Lt21 Lt22 ZERO ; Lt31 Lt32 ZERO]; % inverse of L!  
II = real(L'*Lt');
LLinv= real(Lt'*Lt);


max(max(abs(eye(LSdofs)-II)));

% Testing %

'The maxnorm of the difference between the real part of the generated LU factorization and the stifMatr'
max(max(abs(A-Ah2)));
'Constructing the inverse complex, applying it to the complex total and subtracting Id matr.'
max(max(abs(eye(LSdofs)-(Lt*L))));
spy(Lt*L);

% HOW DOES IT WORK ON A 
%cond(A)
%A = real(Lt2*L*L'*Lt2');

%cond(A)
%A = real(A*Lt');
%cond(A)



%max(abs(A11-real(L11*L11')))
%max(abs(A12-real(L11*L21')))
%max(abs(A13-real(L11*L31')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% In order to approximate the H^-1 norm %%%
%Ah = h^2*Ah;
%fh = h^2*fh;

% Boundary conditions
for I = 1:dofs
	J = 3*I;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N)
    Ah(J,:) = 0;
    Ah(:,J) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
  elseif(j==1 || j==N)
    Ah(J,:) = 0;
    Ah(:,J) = 0;
    Ah(J,J) = 1;
    fh(J) = 0;
  end
end
    
% Boundary conditions
for I = 1:dofs
	J = I+2*dofs;
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  if(i==1 || i==N||j==1 || j==N)
    Ah2(J,:) = 0;
    Ah2(:,J) = 0;
    Ah2(J,J) = 1;
    A(J,:) = 0;
    A(:,J) = 0;
    A(J,J) = 1;
    fh2(J) = 0;
  end
end

uh = Ah\fh;
uh2 = Ah2\fh2;
uh3 = A\fh2;

g = B*uh3;
ERRORUG = norm(U-g(2*dofs+1:end))/norm(U);
if(0)
figure;
surf(x,y,reshape(2*g(1:dofs),N,N));
title('Numerical Solution');
figure;
surf(x,y,reshape(U,N,N));
title('Analytical Solution');

% Plotting
figure;
surf(x,y,reshape(uh2(2*dofs+1:end),N,N));
title('Numerical Solution');

figure;
%surf(x,y,reshape(U,N,N));
%title('Analytical Solution');
surf(x,y,reshape(g(dofs+1:2*dofs),N,N));
title('loading function g part 2 transformed');

figure;
surf(x,y,reshape(fh2(1:dofs),N,N));
title('loading function f part 1 transformed');

figure;
surf(x,y,reshape(g(2*dofs+1:end),N,N));
title('loading function g part 3 transformed');

figure;
surf(x,y,reshape(fh2(dofs+1:2*dofs),N,N));
title('loading function f part 2 transformed');
end

error1 = norm(uh(3:3:LSdofs)-U,'inf')/norm(U,'inf')
error2 = norm(uh2(2*dofs+1:LSdofs)-U,'inf')/norm(U,'inf')
error3 = norm(uh3(2*dofs+1:LSdofs)-U,'inf')/norm(U,'inf')
%conditionNumber1 = cond(Ah)
%conditionNumber2 = cond(Ah2)

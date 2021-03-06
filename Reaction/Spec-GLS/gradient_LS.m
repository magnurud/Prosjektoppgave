function Gh = gradient_LS(mu,B1,B2,W,LDM,dofs)


Bd1 = diag(reshape(B1,dofs,1)); 
Bd2 = diag(reshape(B2,dofs,1)); 

Gh11 = -mu*Bd1*kron(W,W*LDM)  -mu*kron(W,LDM'*W)*Bd1  +Bd1*kron(W,W)*Bd1;
Gh21 = -mu*Bd2*kron(W,W*LDM)  -mu*kron(LDM'*W,W)*Bd1  +Bd1*kron(W,W)*Bd2;
Gh22 = -mu*Bd2*kron(W*LDM,W)  -mu*kron(LDM'*W,W)*Bd2  +Bd2*kron(W,W)*Bd2;
Gh12 = Gh21'; 
ZERO = zeros(size(Gh11));

Gh = sparse([Gh11 Gh12 ZERO; Gh21 Gh22 ZERO ; ZERO ZERO ZERO]);


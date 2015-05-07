function Gh = gradient_LS(mu,B1,B2,W,LDM,dofs)

Gh11 = -mu*B1*kron(W,W*LDM)  -mu*kron(W,LDM'*W)*B1  +B1*kron(W,W)*B1;
Gh21 = -mu*B2*kron(W,W*LDM)  -mu*kron(LDM'*W,W)*B1  +B1*kron(W,W)*B2;
Gh22 = -mu*B2*kron(W*LDM,W)  -mu*kron(LDM'*W,W)*B2  +B2*kron(W,W)*B2;
Gh12 = Gh21'; 
ZERO = zeros(size(Gh11));

Gh = sparse([Gh11 Gh12 ZERO; Gh21 Gh22 ZERO ; ZERO ZERO ZERO]);


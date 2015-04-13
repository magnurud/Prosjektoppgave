function [Ah] = stiffness_2D(N,x,y,wX,wY,LDM);
% description:
%      generate the stiffness matrix Ah;
%
% arguments:
%   - N     Number of points in each direction  
%   - x     x values of the coordinates 
%   - y     y values of the coordinates 
%   - wX    The weight used in the quadrature rule for the x-direction
%   - wY    The weight used in the quadrature rule for the x-direction
%   - LDM   The derivatives of the basisfunctions evaluated in every point
%           ( LDM(i,j) = l_j'(x_i) )
% returns:
%		- Ah 	  Stiffness matrix (dofs^2 elements)
%
% author: Magnus Aa. Rud
% last edit: April 2015

NLS = 3*N; 
dofs = N^2;
LSdofs = NLS^2;
Ah = zeros(LSdofs);
Aloc = zeros(3,3);

for I = 1:dofs
  i = mod(I-1,N)+1;
  j = fix((I-1)/N)+1;
  for J = 1:dofs
    l = mod(J-1,N)+1;
    k = fix((J-1)/N)+1;
		% (1,1)
    if(I==J) 
			Aloc(1,1) = Aloc(1,1) + wX(i)*wY(j);
		end
    if(j==l) 
			for alpha = 1:N
			Aloc(1,1) = Aloc(1,1) + wX(alpha)*wY(j)*LDM(alpha,i)*LDM(alpha,k);
			end
		end
		% (1,2)
		Aloc(1,2) = Aloc(1,2) + wX(j)*wY(k)*LDM(k,i)*LDM(j,l);
		% (1,3)
		if(j==l) 
		Aloc(1,3) = Aloc(1,3) +	wX(i)*wY(j)*LDM(i,k);
		end
		% (2,1)
		Aloc(2,1) = Aloc(2,1) + wX(i)*wY(l)*LDM(l,j)*LDM(i,k);
		% (2,2)
    if(I==J) 
			Aloc(2,2) = Aloc(2,2) + wX(i)*wY(j);
		end
    if(i==k) 
			for beta = 1:N
				Aloc(2,2) = Aloc(2,2) + wX(beta)*wY(j)*LDM(beta,j)*LDM(beta,l);
			end
		end
		% (2,3)
		if(i==k) 
		Aloc(2,3) = Aloc(2,3) +	wX(i)*wY(j)*LDM(j,l);
		end
		% (3,1)
		if(j==l) 
		Aloc(3,1) = Aloc(3,1) +	wX(k)*wY(l)*LDM(k,i);
		end
		% (2,3)
		if(i==k) 
		Aloc(3,2) = Aloc(3,2) +	wX(k)*wY(l)*LDM(l,j);
		end
		% (3,3)
		if(j==l)
			for alpha = 1:N
				Aloc(3,3) = Aloc(3,3) + wX(alpha)*wY(j)*LDM(alpha,i)*LDM(alpha,k);
			end
		end
		if(i==k)
			for beta = 1:N
				Aloc(3,3) = Aloc(3,3) + wX(i)*wY(beta)*LDM(beta,j)*LDM(beta,l);
			end
		end
		Ah((3*I-2):(3*I),(3*J-2):(3*J)) = Ah((3*I-2):(3*I),(3*J-2):(3*J))+Aloc;
	end
end

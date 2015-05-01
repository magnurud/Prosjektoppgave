clc
clear all
close all
format compact
% addpath MWR_Libary

tic
CONE   = 1;%(0/1)
N      = 16;
Nglob  = N*N;
Ymax   = 0.5;
Ymin   = -0.5;
Xmax   = 0.5;
Xmin   = -0.5;
[x,wX] = GLL_(N,Xmin,Xmax);
[y,wY] = GLL_(N,Ymin,Ymax);
T      = 400;
omega  = 2*pi/T;
LagY   = eye(N,N);
LagX   = eye(N,N);
dt     = 1;
dY     = Ymax-Ymin;
dX     = Xmax-Xmin;
THETA  = 0.5;
nROUND = 2;
nSTEP  = T/dt*nROUND;

derY = LagrangeDerivativeMatrix_GLL(N);
derY = derY*2/dY;

derX = LagrangeDerivativeMatrix_GLL(N);
derX = derX*2/dX;


%VELOCITY FIELD
vx = zeros(N);
vy = zeros(N);
[xg,yg] = meshgrid(x,y);
for i=1:N
		for j=1:N
			vx(i,j) = -omega*y(i);
			vy(i,j) = omega*x(j);
		end
end

%PLOT
figure(1)
quiver(vx,vy),xlabel('v_x'),ylabel('v_y')


%INITIAL PROFILE
if CONE==1
		disp('CONE SELECTED')
		fINIT = zeros(N,N);
		for i = 1:N
				for j = 1:N 
						fINIT(i,j) = max( 0, 4*(1-sqrt((i-8)^2+(j-16)^2)/4) );
				end
		end
		fINIT=fINIT';
		Zmax = max(max(fINIT));
else
		disp('BOX SELECTED')
		fINIT = zeros(N,N);
		for i = 1:N
				for j = 1:N 
						if i>5 && i<14 && j>12 && j<20
								fINIT(i,j) = 0.01;
						else
								fINIT(i,j) = 0;
						end
				end
		end
		fINIT=fINIT'; 
		Zmax = max(max(fINIT));
end

%PLOT
figure(2)
subplot(1,2,1),mesh(x,y,fINIT)
xlabel('x'),ylabel('y'),title('Initial profile'),box on
xlim([Xmin,Xmax]),ylim([Ymin,Ymax])


%GATHERING MATRIX
GM = zeros(N,N);
for iY = 1:N
		for iX = 1:N
				iGLOB = iX + N*(iY-1);
				GM(iY,iX) = iGLOB;
		end
end

%GLOBAL WEIGHT
wGLOB = zeros(Nglob,1);
for iY = 1:N
		for iX = 1:N
				iGLOB = GM(iY,iX);
				wGLOB(iGLOB) = wY(iY)*wX(iX);
		end
end
LAMBDA = diag(wGLOB);

f0 = fINIT;

%PROBLEM OPERATOR MATRIX
L = zeros(Nglob,Nglob);
for iY=1:N
		for iX=1:N
				 iGLOB = GM(iY,iX);
				 for jY=1:N
						 for jX=1:N
								 jGLOB = GM(jY,jX);

								 L(iGLOB,jGLOB) = LagX(iX,jX)*LagY(iY,jY) ...
								 -dt*THETA*omega*y(iY)*derX(iX,jX)*LagY(iY,jY) ...
								 +dt*THETA*omega*x(iX)*derY(iY,jY)*LagX(iX,jX);

						 end
				 end
		end
end

A = L'*LAMBDA*L;

%BOUNDAY CONDITIONS 1st STEP - PROBLEM OPERATOR MATRIX
for iY=1:N
		for iX=1:N
				 iGLOB = GM(iY,iX);
				 for jY=1:N
						 for jX=1:N
								 jGLOB = GM(jY,jX);

								 if iY==1 || iY==N || iX==1 || iX==N
										 A(iGLOB,:) = 0;
								 end                 
						 end
				 end
		end
end

%BOUNDARY CONDITIONS 2nd STEP - PROBLEM OPERATOR MATRIX
for iY=1:N
		for iX=1:N
				 iGLOB = GM(iY,iX);
				 for jY=1:N
						 for jX=1:N
								 jGLOB = GM(jY,jX);



								 if iX==1 
										 if vx(iY,iX) >= 0
												A(iGLOB,jGLOB) = LagX(iX,jX)*LagY(iY,jY);
										 else
												A(iGLOB,jGLOB) = derX(iX,jX)*LagY(iY,jY); 
										 end
								 end 


								 if iX==N 
										 if vx(iY,iX) <= 0
												A(iGLOB,jGLOB) = LagX(iX,jX)*LagY(iY,jY);
										 else
												A(iGLOB,jGLOB) = derX(iX,jX)*LagY(iY,jY);
										 end
								 end  


								 if iY==1
										 if vy(iY,iX) >= 0
												A(iGLOB,jGLOB) = LagX(iX,jX)*LagY(iY,jY);
										 else
												A(iGLOB,jGLOB) = derY(iY,jY)*LagX(iX,jX); 
										 end
								 end


								 if iY==N
										 if vy(iY,iX) <= 0 
												A(iGLOB,jGLOB) = LagX(iX,jX)*LagY(iY,jY);
										 else
												A(iGLOB,jGLOB) = derY(iY,jY)*LagX(iX,jX); 
										 end
								end

						 end
				 end
		end
end



%**************************************************************************
%START TIME LOOP
%**************************************************************************


ErrorL2 = zeros(nROUND,1);
PEAK    = zeros(nROUND,1);
fOUT    = zeros(N,N,nROUND);
iROUND  = 1;
disp(['iROUND: ',num2str(iROUND)])
FLAG=0;
k=1;
%=======================
for iSTEP = 1:nSTEP
	%=======================

	if mod(iSTEP,T/dt)==0
		FLAG = 1.0;

		%L2 error
		DUM = sqrt( sum(sum( (fINIT-f0).*(fINIT-f0) )) );
		ErrorL2(k) = DUM;
		disp(['ErrorL2: ',num2str(DUM)])
		%Mass conservation
		DUM = max(max(f0))/max(max(fINIT))*100;
		PEAK(k) = DUM;
		disp(['Conservation of maximum: ',num2str(DUM)])
		%SOLUTION
		fOUT(:,:,k) = f0;
		%UPDATE
		k=k+1;
end
if FLAG==1.0 && mod(iSTEP,T/dt)>0
		iROUND=iROUND+1;
		disp(['iROUND: ',num2str(iROUND)])
		FLAG=0;
end

%SOURCE VECTOR g
g = zeros(Nglob,1);
for iY=1:N
		for iX=1:N
				 iGLOB = GM(iY,iX);

				 dfdx=0;
				 dfdy=0;
				 for jY=1:N
						 for jX=1:N
								 dfdx = dfdx + f0(jY,jX)*derX(iX,jX)*LagY(iY,jY);
								 dfdy = dfdy + f0(jY,jX)*derY(iY,jY)*LagX(iX,jX);
						 end
				 end

				 g(iGLOB) = f0(iY,iX) ...
				 +dt*(1-THETA)*omega*y(iY)*dfdx ...
				 -dt*(1-THETA)*omega*x(iX)*dfdy;


		end
end

F = L'*LAMBDA*g;


%BOUNDARY CONDITIONS - SOURCE VECTOR g
for iY=1:N
		for iX=1:N
				 iGLOB = GM(iY,iX);

				 if iY==1 || iX==1 || iY==N || iX==N
										 F(iGLOB,1) = 0;
								 end 

		end
end

%SOLVE THE SYSTEM MATRIX FORM
fVECTOR = A\F;

%REARRANGE THE SOLUTION VECTOR INTO A MATRIX 
fMATRIC = zeros(N,N);
for iX = 1:N
		for iY = 1:N
				iGLOB = GM(iY,iX);
				fMATRIX(iY,iX) = fVECTOR(iGLOB);

		end
end

%FLOT THE SOLUTION
figure(2)
subplot(1,2,2),mesh(x,y,fMATRIX),colormap([0 0 0])
xlim([Xmin,Xmax]),ylim([Ymin,Ymax]),zlim([0,Zmax])
pause(0.5)


%UPDATE THE SOLTION
f0 = fMATRIX;

end


TIMEused = toc;

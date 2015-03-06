%% Computing the results for 2D poisson equation with Dirichlet Boundary

%% Set up the grid
Nr = 100;
[p, tri, edge] = getDisk(Nr);

%% Create the stiffness matrix and load vector
[Ah] = stiffness_2D(Nr,p,tri);
[fh] = load_2D(Nr,p,tri); %load_2D_reference() also works

%% Imposing the dirichlet boundary
index = length(edge(:,1)) - 1;
for i=(Nr-index):Nr
    fh(i)=0;
    Ah(i,:)=0;
    Ah(:,i)=0; %not nescessary, but convenient (see spy(Ah) with and without)
    Ah(i,i)=1;
end

%% Solving the linear system
uh = Ah\fh;

%% Write the .vtf file for GLview
writeVTF(p,tri,uh,'2DpoissonDirichlet.vtf');

%% Plotting the numerical solution
subplot(1,2,1);
trisurf(tri,p(:,1),p(:,2),uh);
title('Numerical Solution');

%% Plotting the analytical solution
ua=@(x,y) sin(2*pi*(x^2+y^2));
uanal=zeros(Nr,1);
for i=1:Nr
   uanal(i)=ua(p(i,1),p(i,2)); 
end
subplot(1,2,2);
trisurf(tri,p(:,1),p(:,2),uanal);
title('Analytical Solution');

%% Plotting the error
figure(3)
hold on
plot(1:100,uh,'b',1:100,uanal,'r')
xlabel('Node number')
hleg1 = legend('Numerical solution','Analytical solution');
set(hleg1,'Location','NorthEast')
maxerror = max(abs(uh-uanal));
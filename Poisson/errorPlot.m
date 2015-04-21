% errorPlot.m
%
% description:
%      Plotting the error as a funcion of stepsize of the different solutions
%      of the Poisson problem on the square (0,1)^2;
%
% author: Magnus Aa. Rud
% last edit: April 2015

N = [10:1:30];% number of discretization points 
h = ones(1,length(N))./N; %stepsize
eFEM = zeros(1,length(N));
eLSFEM = zeros(1,length(N));

  run FEM/errorFunc.m
  eFEM = ans; 
  run LSFEM/errorFunc.m
  eLSFEM = ans; 

%e = abs((e-e(length(N)))/e(length(N)));
figure;
x = logspace(log(h(end)),log(h(1)));
y = x.^2;
loglog(h,eFEM,'r');
hold on;
loglog(h,eLSFEM,'b');
%hold on;
%loglog(x,y,'b');
legend('FEM','LSFEM')
grid on;


% errorPlot.m
%
% description:
%      Plotting the error as a funcion of stepsize of the different solutions
%      of the Poisson problem on the square (0,1)^2;
%
% author: Magnus Aa. Rud
% last edit: April 2015

N = [10:1:24];% number of discretization points 
h = ones(1,length(N))./N; %stepsize
eFEM = zeros(1,length(N));
eLSFEM = zeros(1,length(N));

  run FEM/errorFunc.m
  eFEM = ans; 
  run LSFEM/errorFunc.m
  eLSFEM = ans; 
  run Spectral/errorFunc.m
  eSpec = ans; 
  run Spectral-LS/errorFunc.m
  eSpecLS = ans; 

%e = abs((e-e(length(N)))/e(length(N)));
figure;
x = logspace(log(h(end)),log(h(1)));
y = x.^2;
loglog(N,eFEM(1,:),'r');
hold on;
loglog(N,eLSFEM(1,:),'b');
legend('FEM','LSFEM')
xlabel('number of nodes in each spacial direction');
ylabel('error');
grid on;
fig1 = gcf;

figure;
loglog(N,eFEM(2,:),'r');
hold on;
loglog(N,eLSFEM(2,:),'b');
legend('FEM','LSFEM')
xlabel('number of nodes in each spacial direction');
ylabel('nondition number');
grid on;
fig2 = gcf;

figure;
loglog(N,eSpec(1,:),'r');
hold on;
loglog(N,eSpecLS(1,:),'b');
legend('Spectral','Spectral LS')
xlabel('Number of nodes in each spacial direction');
ylabel('error');
grid on;
fig3 = gcf;

figure;
loglog(N,eSpec(2,:),'r');
hold on;
loglog(N,eSpecLS(2,:),'b');
legend('Spectral','Spectral LS')
xlabel('number of nodes in each spacial direction');
ylabel('nondition number');
grid on;
fig4 = gcf;


set(fig1, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig1, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig1, '../Latex/Figures/errorFEM-LSFEM', 'pdf') %Save figure

set(fig2, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig2, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig2, '../Latex/Figures/condFEM-LSFEM', 'pdf') %Save figure

set(fig3, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig3, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig3, '../Latex/Figures/errorSpec-SpecLS', 'pdf') %Save figure

set(fig4, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig4, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig4, '../Latex/Figures/condSpec-SpecLS', 'pdf') %Save figure


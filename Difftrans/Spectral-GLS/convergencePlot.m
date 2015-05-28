% ConvergencePlot.m 

N = [10:1:20];% number of discretization points 
h = ones(1,length(N))./N; %stepsize
e1 = zeros(1,length(N));
e2 = zeros(1,length(N));
e3 = zeros(1,length(N));
mu = 1E-2;
alpha = -1;
for i = 1:length(N)
	[e1(i) cn] = runMain_Spec(N(i),mu,alpha); % Spectral
	[e2(i) cn] = runMain_LS(N(i),mu,alpha); % Spectral LS regular
	[e3(i) cn] = runMain_DirFunc(N(i),mu,alpha); %Spectral LS Functional
end

figure;
loglog(N,e1,'r');
hold on;
loglog(N,e2,'b');
loglog(N,e3,'g');
grid on;
xlabel('N')
ylabel('error')
legend('Galerkin','Least-squares: BC Lifting function','Least-squares: BC in functional')

fig1 = gcf;

set(fig1, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig1, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig1, '../../Latex/Figures/Spec_difftrans_Convergence', 'pdf') %Save figure


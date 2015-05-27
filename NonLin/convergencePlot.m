% ConvergencePlot.m 

N = 15;
mu = 1E0;
alpha = 1;
	e1 = runMain_Spec(N,mu,alpha); % Spectral
	e2 = runMain_LS(N,mu,alpha); % Spectral LS regular

figure;
loglog(1:length(e1),e1,'r');
hold on;
loglog(1:length(e2),e2,'b');
grid on;
xlabel('#iterations')
ylabel('error')
legend('Galerkin','Least-squares: BC Lifting function')%,'Least-squares: BC in functional')

fig1 = gcf;

set(fig1, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig1, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig1, '../Latex/Figures/Spec_Nonlin_Convergence', 'pdf') %Save figure



% ConvergencePlot.m 

N = 15;
mu = 1E0;
alpha = -1;
	[e1 cn1 r1] = runMain_Spec_Jnum(N,mu,alpha); % Spectral
	[e2 cn2 r2] = runMain_LS_Jnum(N,mu,alpha); % Spectral LS regular
	[e3 cn3 r3] = runMain_LS_DirFunc_Jnum(N,mu,alpha); % Spectral LS BC-func 

figure;
loglog(1:length(e1),e1,'r');
hold on;
loglog(1:length(e2),e2,'b');
loglog(1:length(e3),e3,'g');
grid on;
xlabel('#iterations')
ylabel('error')
legend('Galerkin','Least-squares: BC Lifting function','Least-squares: BC in functional')

fig1 = gcf;

set(fig1, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig1, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig1, '../Latex/Figures/Spec_Nonlin_Convergence_Jnum', 'pdf') %Save figure

figure;
loglog(1:length(r1),r1,'r');
hold on;
loglog(1:length(r2),r2,'b');
loglog(1:length(r3),r3,'g');
grid on;
xlabel('#iterations')
ylabel('residual')
legend('Galerkin','Least-squares: BC Lifting function','Least-squares: BC in functional')

fig1 = gcf;

set(fig1, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig1, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig1, '../Latex/Figures/Spec_Nonlin_Convergence_Jnum_res', 'pdf') %Save figure




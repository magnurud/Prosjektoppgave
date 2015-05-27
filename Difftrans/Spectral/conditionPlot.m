% ConditionPlot.m 

N = 15;
%alpha = [500:500:5000];% number of discretization points 
alpha = logspace(1,4,20);% number of discretization points 
cn_mu = zeros(1,length(N));
cn_alpha= zeros(1,length(N));
mu = 1;
for i = 1:length(alpha)
	[e cn_alpha(i)] = runMain(N,mu,alpha(i));
	[e cn_mu(i)] = runMain(N,alpha(i),mu);
end

%x = logspace(log(h(1)),log(h(end)));
%x = linspace(h(end),h(1),50);
%y = x.^2;

figure;
%loglog(x,y,'b');
loglog(alpha,cn_alpha,'r');
hold on;
loglog(alpha,cn_mu,'b');
legend('b','mu')
grid on;
fig1 = gcf;
SisteFeil = e(end)
%Kondisjonstall = cn

polyfit(log(alpha),log(cn_mu),1)

set(fig1, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig1, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig1, '../../Latex/Figures/Spec_difftrans_ConditionNumber', 'pdf') %Save figure

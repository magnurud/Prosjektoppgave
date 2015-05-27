% ConvergencePlot.m 

N = [10:1:20];% number of discretization points 
h = ones(1,length(N))./N; %stepsize
e = zeros(1,length(N));
mu = 1E-5;
alpha = -1;
for i = 1:length(N)
	[e(i) cn] = runMain(N(i),mu,alpha);
end

%x = logspace(log(h(1)),log(h(end)));
%x = linspace(h(end),h(1),50);
%y = x.^2;

figure;
%loglog(x,y,'b');
hold on;
loglog(h,e,'r');
xlabel('h')
ylabel('error')
grid on;
fig1 = gcf;
e(end)

%polyfit(h,e,1)

set(fig1, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig1, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig1, '../../Latex/Figures/Spec_difftrans_Convergence', 'pdf') %Save figure

%axis([ax1 ax2],[0 10 -1 1])

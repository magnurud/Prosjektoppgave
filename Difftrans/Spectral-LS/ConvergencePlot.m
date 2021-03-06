% ConvergencePlot.m 

N = [10:1:20];% number of discretization points 
%N = [5:1:15];% number of discretization points 
h = ones(1,length(N))./N; %stepsize
e = zeros(1,length(N));
e2 = zeros(1,length(N));
mu = 1E-5;
alpha = -1;
for i = 1:length(N)
	[e(i) cn] = runMain_DirFunc(N(i),mu,alpha);
	[e2(i) cn] = runMain(N(i),mu,alpha);
end

%x = logspace(log(h(1)),log(h(end)));
%x = linspace(h(end),h(1),50);
%y = x.^2;

figure;
%loglog(x,y,'b');
loglog(h,e,'r');
hold on;
loglog(h,e2,'b');
legend('BC added to functional','BC imposed directly')
xlabel('h')
ylabel('error')
grid on;
fig1 = gcf;
SisteFeil = e(end)
Kondisjonstall = cn

%polyfit(h,e,1)

set(fig1, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(fig1, 'PaperSize', [5.5 5]); %Set the paper to have width 5 and height 5.
saveas(fig1, '../../Latex/Figures/Spec-LS_difftrans_Convergence', 'pdf') %Save figure

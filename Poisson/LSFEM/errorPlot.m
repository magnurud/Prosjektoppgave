% errorPlot.m
%
% description:
%      Plotting the error as a funcion of stepsize of the FEM solution 
%      of the Poisson problem on the square (0,1)^2;
%
% author: Magnus Aa. Rud
% last edit: March 2015

N = [5:5:25];% number of discretization points 
%N = [2 4 8 16 32 64];
h = ones(1,length(N))./N; %stepsize
e = zeros(1,length(N));

for i = 1:length(N)
  n = N(i);
  e(i) = runMain(n);
end

%e = abs((e-e(length(N)))/e(length(N)));

figure;
x = logspace(log(h(end)),log(h(1)));
y = x.^2;
loglog(h,e,'r');
hold on;
loglog(x,y,'b');
legend('loglog of error','reference line with slope = 2')
grid on;


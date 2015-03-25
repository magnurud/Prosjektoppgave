% errorPlot.m
%
% description:
%      Plotting the error as a funcion of stepsize of the FEM solution 
%      of the Poisson problem on the square (0,1)^2;
%
% author: Magnus Aa. Rud
% last edit: March 2015

N = [5:5:25];% number of discretization points 
h = ones(1,length(N))./N; %stepsize
e = zeros(1,length(N));

for i = 1:length(N)
  n = N(i);
  e(i) = runMain(n);
end

%e = abs((e-e(length(N)))/e(length(N)));

loglog(h,e,'r');
refline(2,0);



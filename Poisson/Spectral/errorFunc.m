function eh = errorFunc();
% errorPlot.m
%
% description:
%      Plotting the error as a funcion of stepsize of the FEM solution 
%      of the Poisson problem on the square (0,1)^2;
%
% author: Magnus Aa. Rud
% last edit: March 2015

N = [10:1:24];% number of discretization points 
h = ones(1,length(N))./N; %stepsize
eh = zeros(2,length(N));

for i = 1:length(N)
  n = N(i);
    [eh(1,i) eh(2,i)] = runMain(n);
end


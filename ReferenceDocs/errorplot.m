% This script gives the error plotted against Nr

N = 20;
error = zeros(1,N);
x = zeros(1,N);

for i = 1:N    
    Nr = 500 + 500*(i-1);
    x(i) = Nr;
    error(i) = errorfunction(Nr);
end

%%
plot(x,error,'-x')
xlabel('Number of nodes')
ylabel('Error')
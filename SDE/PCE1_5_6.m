clear all

pdf = @(x,mu,var) 1/sqrt(2*pi*var)*exp(-0.5*(x-mu)^2/(2*var));

M = 10;
N = 20;

tic

ncol = 5*M;
x = linspace(0,M,ncol+1);
y = zeros(1,ncol);
n = 100;

for j = 1:n
    for i = 1:N
        Y(i) = M*rand();
    end
    A(j) = 1/N*sum(Y);
    for i = 1:ncol+1
        if(A(j)>=x(i) && A(j)<x(i+1))
            y(i) += 1;
        end
    end
end

toc

bar(x(1:end-1),y);

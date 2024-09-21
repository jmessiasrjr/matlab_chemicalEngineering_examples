clear all

n = 10000;
nsubi = 20;

tic

s = linspace(0,1,nsubi+1);
y = zeros(1,nsubi);

for i = 1:n
    x(i) = rand();
    for j = 1:nsubi
        if(x(i)>s(j) && x(i)<=s(j+1))
            y(j) += 1;
        end
    end
end

mean = 1/n*sum(x);
for i = 1:n
    X(i) = (x(i)-mean)^2;
end
var = 1/(n-1)*sum(X);

toc

bar(y);

fprintf('mean    : %-10.6f\n',mean);
fprintf('variance: %-10.6f\n',var);
        

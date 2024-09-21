clear all

n = 5000;
nsubi = 50;

lamb = 0.5;
N = 5;
a = 2;

tic

s = linspace(0,N,nsubi+1);
for j = 1:nsubi
    Xb(j) = (s(j+1)+s(j))/2;
end
y = zeros(1,nsubi);
k = 1;

for i = 1:n
    x(i) = -log(1 - rand())/lamb;
    for j = 1:nsubi
        if(x(i)>s(j) && x(i)<=s(j+1))
            y(j) += 1;
        end
    end
    if(x(i)>=a)
        xx(k) = x(i);
        k += 1;
    end
end

mean = 1/n*sum(x);
for i = 1:n
    X(i) = (x(i)-mean)^2;
end
var = 1/(n-1)*sum(X);

toc

bar(Xb,y);
xlabel('x');ylabel('N');

fprintf('lambda = %-10.6f\n',lamb);
fprintf('mean    : %-10.6f\n',mean);
fprintf('variance: %-10.6f\n',var);
fprintf('mean X|x>=%-10.6f = %-10.6f\n',a,sum(xx/length(xx)));

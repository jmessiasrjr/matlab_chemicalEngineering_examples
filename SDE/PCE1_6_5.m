clear all

tic

a = 0.1;
b = 0.01;

Prob = @(a,b) [1-a a;b 1-b];

ny = 50;
N = 200;
P = Prob(a,b);
city = zeros(N,ny+1);

cs = [0 1];% E=0; W=1

city(:,1) = cs(1);
tE = 0;

for j = 1:N
    for i = 2:(ny+1)
        p = rand();
        if(city(j,i-1)==cs(1) && p>=P(1,1))
            city(j,i)=cs(2);
        elseif(city(j,i-1)==cs(2) && p>=P(2,2))
            city(j,i)=cs(1);
        else
            city(j,i)=city(j,i-1);
        end
    end
    if(city(j,end)==cs(1)) tE +=1; end
end

sc = sum(city);
avE = sc;
avW = abs(N-sc);

toc

plot((0:ny),avE,'linewidth',2,'-k');hold on
plot((0:ny),avW,'linewidth',2,'-r');hold off
title(sprintf('Persons living in city %c and %c\n',char(69),char(87)));
xlabel('years');ylabel('persons');
legend('E','W','location','east');
ylim([0 N]);

fprintf('Persons living in E = %d\n',tE);
fprintf('Persons living in W = %d\n',N-tE);

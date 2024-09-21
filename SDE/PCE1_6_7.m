clear all

tic

a = 0.1;
b = 0.01;
T = 500;
n = 100;

Prob = @(a,b) [1-a a;b 1-b];

p = 1./11;
P = Prob(a,b);

fX = [1 -1];% [E ; W]
pp = [p 1-p];

freq = zeros(n,T+1);
%values = zeros(n,T+1);
X = zeros(1,T+1);

for i = 1:n
    X(1) = GRNtp(pp(1),fX(1),fX(2));
    E = zeros(1,T+1);
    if(X(1)==fX(1)) E(1)=1; end
        for t = 1:T
            if(X(t)==fX(1))
                X(t+1)=GRNtp(P(1,1),fX(1),fX(2));
            else
                X(t+1)=GRNtp(P(2,1),fX(1),fX(2));
            end
            if(X(t+1)==fX(1)) E(t+1)=1; end
        end
    freq(i,:) = E;
    %values(i,:) = X;
end

nE = sum(freq);
nW = n*ones(1,T+1) - nE;

toc

figure(1);
plot((0:T),nE,'linewidth',2,'-k');hold on
plot((0:T),nW,'linewidth',2,'-r');hold off
title(sprintf('Population distribution of %c and %c\n',char(69),char(87)));
xlabel('weeks');ylabel('persons');
legend('E','W','location','east');

fprintf('Relative frequency of persons living in E after %-d weeks = %-10.4f\n',T,nE(end)/n);
%figure(2);surf((0:T),(1:n),values);xlabel('weeks');ylabel('person');zlabel('X_i'); %If there are more than 20 persons or 100 time steps, too much memory

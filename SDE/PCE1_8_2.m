clear all

tic

p = 0.5;

T0 = 0; T = 1.;
X0 = 0.;

n = 100;
m = 2; % join path

x = [-1 1];
DT = (T-T0)/n;

SQDT = sqrt(DT);
TK = zeros(1,n+1);
SNT = zeros(1,n+1);
% IC
TK(1) = T0; 
TK2(1) = T0; 
SK = 0;SK2 = 0;
SNT(1) = 0;
SNT2(1) = 0;
k = 1;

for i = 2:n+1
    TK(i) = TK(i-1) + DT;
    X(i-1) = GRNtp(p,x(1),x(2));
    SK = SK + X(i-1);
    SNT(i) = SK*SQDT;
    if(mod(i,m)==m-1)
        SK2 = SK2 + sum(X(max(1,i-m):i-1));
        TK2(k+1) = TK2(k) + min(i-1,m)*DT; 
        SNT2(k+1) = SK2*sqrt(min(i-1,m))*SQDT;
        k += 1;
    end
end

toc

plot(TK,SNT,'linewidth',2,'-k');hold on;
plot(TK2,SNT2,'linewidth',4,'-k');hold off;
grid;
xlabel('time');ylabel('W(t)');

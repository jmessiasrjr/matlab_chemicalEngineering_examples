clear all

tic

p = 0.5;

T0 = 0; T = 1.;
X0 = 0.;

n = 100;

x = [-1 1];
DT = (T-T0)/n;

SQDT = sqrt(DT);
TK = zeros(1,n+1);
SNT = zeros(1,n+1);
% IC
TK(1) = T0; 
SK = 0;
SNT(1) = 0;

for i = 2:n+1
    TK(i) = TK(i-1) + DT;
    X = GRNtp(p,x(1),x(2));
    SK = SK + X;
    SNT(i) = SK*SQDT;
end

toc

plot(TK,SNT,'linewidth',2,'-k');
xlabel('time');ylabel('W(t)');

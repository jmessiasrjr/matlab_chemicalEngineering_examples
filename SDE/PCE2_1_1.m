clear all

tic

n = 2000;
T0 = 0; T = 5;
DT = (T - T0)/n;

SQDT = sqrt(DT);
TK = zeros(1,n+1);
WT = zeros(1,n+1);
XT = zeros(1,n+1);
RST = zeros(1,n+1);

TK = T0;WT = 0.;XT = 0.;RST = 0.;

for i = 2:n+1
    TK(i) = TK(i-1) + DT;
    if(mod(i,2)==0)
        [X1 X2] = GRNnorm();
    else
        X1 = X2;
    end
    DWT = X1*SQDT;
    WT(i) = WT(i-1)+DWT;
    XT(i) = 0.5*(WT(i)*WT(i)-TK(i));
    RST(i) = RST(i-1)+(WT(i)-DWT)*DWT;
end

toc

plot(TK,XT,'linewidth',2,'-k');hold on;
plot(TK,WT,'linewidth',2,':b');
plot(TK,RST,'linewidth',2,'-r');hold off;
legend('Exact Ito','Wiener process','Ito sum');
xlabel('time');ylabel('X(t)');

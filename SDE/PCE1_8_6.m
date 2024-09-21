clear all

tic

p = 0.5;

T0 = 0; T = 1.;
X0 = 0.;
x = [-1 1];

n = 100;

% IC
TK(1) = T0; 
SK = 0;SK2 = 0;
SNT(1) = 0;

num = round(n/2+.1);
DT = (T-T0)/n;

SQDT = sqrt(DT);
TK = zeros(1,n+1);
SNT = zeros(1,n+1);
ratio = zeros(1,num);

k = 1;
h = 1.e-7;
TA = h;

for i = 2:n+1
    TK(i) = TK(i-1) + DT;
    X(i-1) = GRNtp(p,x(1),x(2));
    SK = SK + X(i-1);
    SNT(i) = SK*SQDT;
    if(i>=num+1)
        ratio(k+1) = (SNT(i)-SNT(i-1))/h;
        h += DT;
        TA(k+1) = h;
        k += 1;
    end
end

toc

plot(TK,SNT,'linewidth',2,'-k');hold on;
plot((T-T0)/2+TA,ratio,'linewidth',2,'-r');hold off;
grid;
xlabel('time');ylabel('W(t)');
ylim([-1.1*max(abs(SNT)) 1.1*max(abs(SNT))]);
legend('W(t)','W(t)/t');

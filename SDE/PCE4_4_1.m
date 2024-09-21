clear all

tic

T0 = 0; T = 1.;
X0 = 1.; WPT0 = 0.;

a = 1; b = 2.;

n = 2000;

DT = (T-T0)/n;

SQDT = sqrt(DT);
XT1 = X0*ones(1,n+1);
XT2 = X0*ones(1,n+1);
TK = zeros(1,n+1);
% IC
TK(1) = T0; WT = WPT0;

for i = 2:n+1
    TK(i) = TK(i-1) + DT;
    if(mod(i,2)==0)
        [X1 X2] = GRNnorm();
    else
        X1 = X2;
    end
    WT = WT + X1*SQDT;
    XT1(i) = X0*exp( (a-0.5*b*b)*(TK(i)-T0) + b*(WT-WPT0) );
    XT2(i) = (WT + sqrt(X0))^2;
end

toc

plot(TK,XT1,'linewidth',2,'-k');hold on;
plot(TK,XT2,'linewidth',2,'-r');hold off;
xlabel('time');ylabel('X(t)');

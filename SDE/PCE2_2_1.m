clear all

tic

T0 = 0; T = 1.;
X0 = 1.; WPT0 = 0.;

a = 1; b = 1.;

n = 2000;

DT = (T-T0)/n;

SQDT = sqrt(DT);
XT = ones(1,n+1);
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
    XT(i) = X0*exp( (a-0.5*b*b)*(TK(i)-T0) + b*(WT-WPT0) );
end

toc

plot(TK,XT,'linewidth',2,'-k');
xlabel('time');ylabel('X(t)');

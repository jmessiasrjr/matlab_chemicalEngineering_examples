clear all

tic

T0 = 0; T = 1.;
X0 = 0.;

n = 1e4;

DT = (T-T0)/n;

SQDT = sqrt(DT);
WT = zeros(1,n+1);
TK = zeros(1,n+1);
% IC
TK(1) = T0; WT(1) = 0;

for i = 2:n+1
    TK(i) = TK(i-1) + DT;
    if(mod(i,2)==0)
        [X1 X2] = GRNnorm();
    else
        X1 = X2;
    end
    WT(i) = WT(i-1) + X1*SQDT;
end

toc

plot(TK,WT,'linewidth',2,'-k');
xlabel('time');ylabel('W(t)');

EWX = sum(WT)/n;

fprintf('E(W(X)) = %-10.6f\n',EWX);

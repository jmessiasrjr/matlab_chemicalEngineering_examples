clear all

n = 10000;
h = 10;

tic

Cov = [h h*h/2;h*h/2 h*h*h/3];
EX1 = 0; EX2 = 0.;
E2X1 = 0.; E2X2 = 0;
EX1X2 = 0.;

for i = 1:n
    [G1 G2] = GRNcorr(Cov);
    EX1 += G1; EX2 += G2;
    E2X1 += G1*G1; E2X2 += G2*G2;
    EX1X2 += G1*G2;
end

EX1 = EX1/n; EX2 = EX2/n;
E2X1 = E2X1/n; E2X2 = E2X2/n;
EX1X2 = EX1X2/n;

toc

fprintf('EX1   : %-10.6f E(X1)   : %-10.6f\n',EX1,0.);
fprintf('EX2   : %-10.6f E(X2)   : %-10.6f\n',EX2,0.);
fprintf('E2X1  : %-10.6f E(X1^2) : %-10.6f\n',E2X1,Cov(1,1));
fprintf('E2X2  : %-10.6f E(X2^2) : %-10.6f\n',E2X2,Cov(2,2));
fprintf('EX1X2 : %-10.6f E(X1X2) : %-10.6f\n',EX1X2,Cov(1,2));

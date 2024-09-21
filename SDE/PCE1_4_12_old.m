clear all

fdp = @(detC,summ) sqrt(detC)/(2*pi)*exp(-1/2*(summ));

% RANDU generator IBM
a = 16807; b = 0; c = 2^31-1;

tic

n = 2;
h = 1;

mu = [0 0];

C = [h h^2/2;h^2/2 h^3/3];
detC = det(C);
iC = inv(C);

for i = 1 : n
    ni = randi([0,c-1]);
    x = mod(a*ni+b,c);
    U1 = x/c; % RANDU generator IBM
    U2 = rand(); % MATLAB random generator number
    X(1,i) = sqrt(-2*log(U1))*cos(2*pi*U2);% Box-Muller
    X(2,i) = sqrt(-2*log(U1))*sin(2*pi*U2);
    Xm = [X(1,i)-mu(1);X(2,i)-mu(2)];
    cij = Xm'*C*Xm;
    p(i) = fdp(detC,cij);
end

toc

plot3(X(1,:),X(2,:),p,'linewidth',2,'ok');
xlabel('X1');ylabel('X2');zlabel('p');

mean1 = 1/n*sum(X(1,:));
mean2 = 1/n*sum(X(2,:));

for i = 1:n
    X12(i) = X(1,i)*X(2,i);
    V1(i) = (X(1,i)-mean1)^2;
    V2(i) = (X(2,i)-mean2)^2;
    V12(i) = ((X(1,i)+X(2,i))-(mean1+mean2))^2;
end

EX1X2 = 1/n*sum(X12);

var1 = 1/(n-1)*sum(V1);
var2 = 1/(n-1)*sum(V2);
var12 = 1/(n-1)*sum(V12);

fprintf('mean of X1    : %-10.6f\n',mean1);
fprintf('mean of X2    : %-10.6f\n',mean2);
fprintf('Var(X1)       : %-10.6f\n',var1);
fprintf('Var(X2)       : %-10.6f\n',var2);
fprintf('E(X1*X2)      : %-10.6f | %10.6f\n',EX1X2,mean1*mean2);
fprintf('Var(X1 + X2)  : %-10.6f | %10.6f\n',var12,var1+var2);

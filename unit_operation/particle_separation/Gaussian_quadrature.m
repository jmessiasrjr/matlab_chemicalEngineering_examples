clear all

% Gauss-Legendre quadrature
    x = [-0.973906 -0.865063 -0.679409 -0.433395 -0.148874 0.148874 0.433395 0.679409 0.865063 0.973906];
    w = [0.066671 0.149451 0.219086 0.269267 0.295524 0.295524 0.269267 0.219086 0.149451 0.066671];

Q = 27.7/60;
n = 1.5;
K = 17.3;
k = 0.095;
rhos = 2500;
rhof = 1.1;
mu = 1.7e-5;
Dc1 = 0.636;
Dc2 = 0.45;
Cv = 0.03;
eps = 1 - Cv;

U = @(Dc,E) 8*Q/(Dc*Dc)*(2*E-1);
DS = @(Dc, U, E) 1e6*Dc*k*sqrt(mu/(0.125*U*Dc*(2*E-1)*(rhos-rhof)));
DD = @(x) K*(log(1/(1-x)))^(1/n);
G = @(DDstar) (DDstar)^2/(1 + (DDstar)^2);

U1 = U(Dc1,eps);
U2 = U(Dc2,1.);
Ds1 = DS(Dc1,U1,eps);
Ds2 = DS(Dc2,U2,1.);

for i = 1:10
    f(i) = DD(0.5*x(i)+0.5);
    g1(i) = G(f(i)/Ds1);
    F1(i) = w(i)*g1(i);
end
Et1 = 0.5*sum(F1);
for i = 1:10
    f(i) = DD(0.5*x(i)+0.5);
    g2(i) = (1. - G(f(i)/Ds1))*G(f(i)/Ds2);
    F2(i) = w(i)*g2(i);
end
Et2 = 0.5*sum(F2);

fprintf('Cut diameter 1: \t%-8.4f µm\t Et: %-10.4f\n',Ds1,Et1);
fprintf('Cut diameter 2: \t%-8.4f µm\t Et: %-10.4f\n',Ds2,Et2);
fprintf('Total efficiency: %-10.4f\n',Et1+Et2);

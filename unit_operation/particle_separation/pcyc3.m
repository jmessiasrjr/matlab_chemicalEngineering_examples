clear all

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
X = @(D) 1 - exp(-(D/K)^n);
U = @(Dc,E) 8*Q/(Dc*Dc)*(2*E-1);
DS = @(Dc, U, E) 1e6*Dc*k*sqrt(mu/(0.125*U*Dc*(2*E-1)*(rhos-rhof)));
I = @(Ds) ((1.11*n/(0.118+n))/(1.81-0.322*n+(K/Ds)))*(K/Ds);

U1 = U(Dc1,eps)
U2 = U(Dc2,1.)
Ds1 = DS(Dc1,U1,eps)
Ds2 = DS(Dc2,U2,1.)
Et1 = I(Ds1)
Et2 = (1-Et1)*I(Ds2)

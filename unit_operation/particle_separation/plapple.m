clear all

Q = 100./60;
n = 1.5;
K = 37.7;
k = 0.095;
rhos = 2300;
rhof = 4.43;
mu = 3.5e-5;
Et = 0.85;
u = 15;

X = @(D) 1 - exp(-(D/K)^n);
I = @(Ds) ((1.11*n/(0.118+n))/(1.81-0.322*n+(K/Ds)))*(K/Ds);
P = @(Ds,Dc) k*sqrt(mu/(0.125*u*Dc*(rhos-rhof)))-Ds/Dc;
P2 = @(Ds,Dc,Q) k*sqrt(mu*Dc/(Q*(rhos-rhof)))-Ds/Dc;

Dstar = 1e-6*fzero(@(Ds) I(Ds)-Et,1)

Dcyl = fzero(@(Dc) P(Dstar,Dc),0.01)
Q1 = u*Dcyl*Dcyl/8;

N = Q/Q1
Q2 = Q/ceil(N);
Dcyln = fzero(@(Dc) P2(Dstar,Dc,Q2),0.01)
unew = 8*Q2/(Dcyln*Dcyln)

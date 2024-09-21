clear all

Dp = [20 30 40 50 60 70 80 100];
X = [.15 .28 .48 .54 .64 .72 .78 .88];

x = polyfit(Dp,X,1);

b = 9.81;
rhof = 1000;
mu = 1.e-3;
rhos1 = 2200;
rhos2 = 2600;
phi1 = 0.6;
phi2 = 0.8;
B = 3; H = 0.3; L = 4;

k1 = @(phi) 0.843*log10(phi/0.065);
k2 = @(phi) 5.31-4.88*phi;
K11 = k1(phi1);
K21 = k2(phi1);
K12 = k1(phi2);
K22 = k2(phi2);
VT = @(D, rs, K1) (rs - rhof)*b*D*D*K1/(18*mu);

vts = VT(70e-6,rhos2,K12);
um = L/H*vts;
Q = um*H*B;

Re = @(Dp, U) Dp*U*rhof/mu;
CdpRe = @(U, rs) 4./3*(rs-rhof)*mu*b/(rhof^2*(U^3));     %Cd/Re
fRe = @(K1,K2,CdpRe) ((24./(K1*CdpRe))^(1.3/2) + (K2/CdpRe)^(1.3))^(1./1.3);

cdpre = CdpRe(vts, rhos1);
re1 = fRe(K11,K12,cdpre);
dp1 = fzero(@(Dp) Re(Dp, vts)-re1,0);
xf = polyval(x,dp1*1e6);

fprintf('Volumetric flow: \t\t\t\t\t%-10.6f m^3/s\n',Q*3600);
fprintf('Diameter of the larger particle of chalk on sand: \t%-10.6f µm\n',dp1*1e6);
fprintf('Lost Percentage: \t\t\t\t\t%-10.6f %%\n',(1.-xf)*100);

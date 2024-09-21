clear all

Dt = [.03 .04 .06 .12];
m = [0.00462 0.00675 0.00775 0.00442];
b = 9.81;

mt = 0.025;
Q = 37e-6/60;
rhof = 1000;
mu = 1e-3;
rhos = 1800;

vel = @(Dt) 4*Q/(pi*Dt*Dt);

K1 = @(phi) 0.843*log10(phi/0.065);
K2 = @(phi) 5.31-4.88*phi;

CdpRe = @(U) 4./3*(rhos-rhof)*mu*b/((rhof^2)*(U^3));     %Cd/Re
Re = @(K1,K2,CdpRe) ((24./(K1*CdpRe))^(1.3/2) + (K2/CdpRe)^(1.3))^(1./1.3);
dp = @(re, U) re*mu/(U*rhof)*1e6;
X = @(M) (mt-M)/mt;

for i = 1:4
    v(i) = vel(Dt(i));
    cdpre(i) = CdpRe(v(i));
    re(i) = Re(1,0.43,cdpre(i));
    DP(i) = dp(re(i),v(i));
    x(i) = X(sum(m(1:i)));
    fprintf('Elutriador %-d: v [m/s] = %-10.6f Dp [µm] = %-10.6f X = %-10.6f\n',i,v(i),DP(i),x(i));
end


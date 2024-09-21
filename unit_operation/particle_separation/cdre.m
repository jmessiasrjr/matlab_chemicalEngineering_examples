clear all

tic

rhosa = 2200; phia = 0.70;
rhosb = 3200; phib = 0.85;
rhof = 1000; mu = 1e-3;
b = 9.81;

Dpmin = 0.149e-3; Dpmax = 0.595e-3;

Re = @(Dp, U) Dp*U*rhof/mu;
CdRe2 = @(Dp, rhos) 4./3*(rhof*(rhos-rhof)*(Dp^3)*b)/(mu*mu); %Cd*Re^2
CdpRe = @(U, rhos) 4./3*(rhos-rhof)*mu*b/(rhof^2*(U^3));     %Cd/Re

K1 = @(phi) 0.843*log10(phi/0.065);
K2 = @(phi) 5.31-4.88*phi;

% Coelho and Massarani correlation 0.65<phi<=1 Re<5e4
fCd = @(K1,K2,Re) ((24./(K1*Re))^(0.85) + K2^(0.85))^(1./0.85);
fRe = @(K1,K2,CdRe2) ((24./(K1*CdRe2))^(1.2) + (K2/CdRe2)^(1.2/2))^(-1./1.2);
fRe2 = @(K1,K2,CdpRe) ((24./(K1*CdpRe))^(1.3/2) + (K2/CdpRe)^(1.3))^(1./1.3);

Dp = Dpmin;
i = 1;
while (Dp < Dpmax)
    x(i) = Dp;
    k1a = K1(phia); k1b = K1(phib);
    k2a = K2(phia); k2b = K2(phib);
    cdre2a = CdRe2(Dp,rhosa);
    cdre2b = CdRe2(Dp,rhosb);
    rea = fRe(k1a,k2a,cdre2a);
    reb = fRe(k1b,k2b,cdre2b);
    ua(i) = fzero(@(U) Re(Dp, U)-rea,0);
    ub(i) = fzero(@(U) Re(Dp, U)-reb,0);
    i = i + 1;
    Dp = Dp + 1e-6;
end

toc

plot(x,ua,'linewidth',2,'-k');
hold on;
plot(x,ub,'linewidth',2,'-b');
hold off;
legend('leve','pesado');
xlabel('Particle diameter [m]');
ylabel('Terminal velocity [m/s]');
grid;

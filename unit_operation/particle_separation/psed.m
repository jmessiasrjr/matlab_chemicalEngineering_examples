clear all

function uvinf = CA(eps, Reinf)
    if(Reinf < 0.2)
        if(eps > 0.9)
            uvinf = 4.8*eps - 3.8;
        else
            uvinf = 0.83*eps^(3.94);
        end
    elseif(Reinf < 500)
        A = 0.28*eps^(-5.96);
        B = 0.35 - 0.33*eps;
        uvinf = 1./(1.+A*Reinf^(-B));
    else
        uvinf = 0.095*exp(2.29*eps);
    end
end

b = 9.81;
rhof = 1.2;
mu = 1.8e-5;
rhos = 2500;
phi = 0.7;
B = 2; H = 2; L = 16; d = 0.1;
Q = 4; %[m^3/s]
Cv1 = 0.002;
Cv2 = 0.05;

eps = @(Cv) 1 - Cv;
k1 = @(phi) 0.843*log10(phi/0.065);
k2 = @(phi) 5.31-4.88*phi;
K1 = k1(phi);
K2 = k2(phi);
VT = @(D) (rhos - rhof)*g*D*D*K1/(18*mu);

um = Q/(H*B);
vt = d/L*um;

Re = @(Dp, U) Dp*U*rhof/mu;
CdpRe = @(U) 4./3*(rhos-rhof)*mu*b/(rhof^2*(U^3));     %Cd/Re
fRe = @(CdpRe) ((24./(K1*CdpRe))^(1.3/2) + (K2/CdpRe)^(1.3))^(1./1.3);

vinf = @(cv) vt*(1. + 2.5*cv);

vinf1 = vinf(Cv1);
cdpre1 = CdpRe(vinf1);
re1 = fRe(cdpre1);
dp1 = fzero(@(Dp) Re(Dp, vinf1)-re1,0);

vinf2 = vinf(Cv2);
cdpre2 = CdpRe(vinf2);
re2 = fRe(cdpre2);
dp2 = fzero(@(Dp) Re(Dp, vinf2)-re2,0);

fprintf('Terminal velocity: %-8.6f m\n',vt);
fprintf('Particle diameter: %-10.8f µm\n',dp1*1e6);
fprintf('Particle diameter: %-10.8f µm\n',dp2*1e6);

clear all

% Richardson-Zaki
function uvinf = RZ(eps, Reinf)
    if(Reinf < 0.2)
        uvinf = eps^(3.65);
    elseif(Reinf < 1)
        uvinf = eps^(4.35*Reinf^(-0.03)-1);
    elseif(Reinf < 500)
        uvinf = eps^(4.45*Reinf^(-0.1)-1);
    else
        uvinf = eps^(1.39);
    end
end

% Concha-Almendra
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

rhof = .9; mu = 2.3e-2;
rhos = 2.3; dp = 0.08; phi = .8;
Cm = 0.26;
g = 981;

K1 = @(phi) 0.843*log10(phi/0.065);
K2 = @(phi) 5.31-4.88*phi;

eps = 1 - Cm/rhos; %porosity
CdRe2 = 4./3*(rhof*(rhos-rhof)*(dp^3)*g)/(mu*mu);
Reinf = ((24./(K1(phi)*CdRe2))^(1.2) + (K2(phi)/CdRe2)^(1.2/2))^(-1./1.2);
vinf = mu*Reinf/(rhof*dp); %v_{\infty}
uvinf = CA(eps, Reinf); 
U = uvinf*vinf; %||u-v||

fprintf('porosity: %-6.4f  v_{inf}: %-8.6f  U: %-8.6f\n',eps,vinf,U);

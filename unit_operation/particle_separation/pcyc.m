clear all

b = 9.81;
K = 37.7;
k = 0.095;
n = 1.5;
rhof = 1.028;
mu = 2.052e-5;
rhos = 1050;
Dc = 0.55;
Bc = 0.14/Dc;
Hc = 0.275/Dc;
u = 15;

DS = 1e6*k*Dc*sqrt(mu/((Bc*Hc*Dc*u)*(rhos-rhof)));
eta = @(D) ((D/DS)^2)/(1+((D/DS)^2));

for i = 1:500
    Ds(i) = 1.e-1*(i);
    Ei(i) = eta(Ds(i));
end

fprintf('Efficiency for 20 µm: %-10.6f\n',eta(20));
plot(Ds,Ei,'linewidth',2,'-k');
xlabel('Diameter [10^{-6} m]');
ylabel('Individual Efficiency');grid;

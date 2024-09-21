clear all

b = 9.81;
k = 0.095;
N = 3;
rhof = 0.7459;
mu = 2.577e-5;
rhos = 3000;
Dc = 0.508;
Bc = 0.25;
Hc = 0.5;
u = 15;
beta = 315;

Q = 3*u*Hc*Bc*Dc*Dc;
DS = 1e6*k*Dc*sqrt(mu/((Bc*Hc*Dc*u)*(rhos-rhof)));
eta = @(D) ((D/DS)^2)/(1+((D/DS)^2));
Uc = 4*Q/(3*pi*Dc*Dc);
deltap = @(U) beta*rhof*U*U/2;
P = @(delp) delp*Q/(746*0.5);

for i = 1:500
    Ds(i) = 1.e-1*(i);
    Ei(i) = eta(Ds(i));
end

Dp = fzero(@(D) eta(D)-0.95,[0 50]);
delp = deltap(Uc);
pot = P(delp);

fprintf('Capacity: \t\t\t\t\t%-10.6f m^3/s\n',Q);
fprintf('Particle diameter with efficiency of 95%%: \t%-10.6f µm\n',Dp);
fprintf('Potency (50%% efficiency): \t\t\t%-10.6f HP\n',pot);
plot(Ds,Ei,'linewidth',2,'-k');
xlabel('Diameter [10^{-6} m]');
ylabel('Individual Efficiency');grid;

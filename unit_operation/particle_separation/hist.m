clear all

rhos = 3000;
rhofar = 1.2; rhofh2o = 1000;
muar = 1.8e-5; muh2o = 1.e-3;
Dp1 = 43e-6; Dp2 = 77e-6;

A = @(rhof,mu,Dp) 36*mu/((2*rhos + rhof)*Dp^2);
vpvt = @(A,t) (exp(A*t) - 1)/(exp(A*t));

A1 = A(rhofar,muar,Dp1);
A2 = A(rhofh2o,muh2o,Dp2);

for i = 1:501
    t(i) = (i-1)/10000;
    vpvt1(i) = vpvt(A1,t(i));
    vpvt2(i) = vpvt(A2,t(i));
end

plot(t,vpvt1,'linewidth',2,'-k',t,vpvt2,'linewidth',2,'-r');
legend('ar','agua','location','southeast');
xlabel('t [s]'); ylabel('v/v_t'); grid;

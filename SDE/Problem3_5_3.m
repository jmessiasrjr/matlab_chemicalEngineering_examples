clear all

a1 = 1; a2 = 15;
a = @(alpha,x,t) -alpha*x;
Xe1 = @(x0,t) x0*exp(-a1*t);
Xe2 = @(x0,t) x0*exp(-a2*t);

tic

n = 128;
t0 = 0.; T = 1.;
x0 = 1.;

te = linspace(0,T-t0,1000);
X1 = Xe1(x0,te);
X2 = Xe2(x0,te);

dt = (T-t0)/n;
x1 = x0*ones(1,n+1);
x2 = x0*ones(1,n+1);
t = t0*ones(1,n+1);

for i = 2:n+1
    t(i) = t(i-1) + dt;
    x1(i) = x1(i-1) + a(a1,x1(i-1),t(i-1))*dt;
    x2(i) = x2(i-1) + a(a2,x2(i-1),t(i-1))*dt;
end

toc

figure(1)
plot(te,X1,'linewidth',4,'-k');hold on;
plot(te,X2,'linewidth',4,'color',[0.2 0.2 0.2]);
plot(t,x1,'linewidth',2,'-r');
plot(t,x2,'linewidth',2,'-b');hold off;
legend('alpha1','alpha2','exact1','exact2');
xlabel('time');ylabel('x(t)');
figure(2)
plot3(t,x1,x2,'linewidth',2,'-k');
xlabel('time');ylabel('real');zlabel('imaginary');

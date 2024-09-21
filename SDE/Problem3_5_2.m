clear all

alpha = 1;
a = @(x,t) -alpha*x;
Xe = @(x0,t) x0*exp(-alpha*t);

tic

n = 8;
t0 = 0.; T = 1.;
x0 = 1.;

te = linspace(0,T-t0,1000);
X = Xe(x0,te);

dt = (T-t0)/n;
xh = x0*ones(1,n+1);
xt = x0*ones(1,n+1);
t = t0*ones(1,n+1);

for i = 2:n+1
    t(i) = t(i-1) + dt;
    xh(i) = xh(i-1) + 0.5*( a(xh(i-1),t(i-1)) + ...
            a(xh(i-1) + a(xh(i-1),t(i-1))*dt) )*dt;
    xt(i) = (1 - 0.5*alpha*dt)/(1 + 0.5*alpha*dt)*xt(i-1);
end

toc

plot(t,xh,'linewidth',2,'-r');hold on;
plot(t,xt,'linewidth',2,'-b');
plot(te,X,'linewidth',4,'-k');hold off;
legend('Heun','trapezoidal','exact');
xlabel('time');ylabel('x(t)');

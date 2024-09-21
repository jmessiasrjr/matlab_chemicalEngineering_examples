clear all

a = @(x,t) -5*x;
Xe = @(x0,t) x0*exp(-5*t);

tic

n = 8;
t0 = 0.; T = 1.;
x0 = 1.;

te = linspace(0,T-t0,1000);
X = Xe(x0,te);

dt = (T-t0)/n;
x = x0*ones(1,n+1);
t = t0*ones(1,n+1);

for i = 2:n+1
    x(i) = x(i-1) + a(x(i-1),t(i-1))*dt;
    t(i) = t(i-1) + dt;
end

toc

plot(t,x,'linewidth',2,'-r');hold on;
plot(te,X,'linewidth',2,'-k');hold off;
legend('numerical','exact');
xlabel('time');ylabel('x(t)');

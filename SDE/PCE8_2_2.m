clear all

a = @(x,t) t*x*(2 - x);

tic

n = 8;
t0 = 0.; T = .5;
x0 = 1.;

dt = (T-t0)/n;
x = x0*ones(1,n+1);
xrk = x0*ones(1,n+1);
t = t0*ones(1,n+1);

for i = 2:n+1
    t(i) = t(i-1) + dt;
    x(i) = x(i-1) + a(x(i-1), t(i-1))*dt;
    kn1 = a(xrk(i-1), t(i-1));
    kn2 = a(xrk(i-1) + 0.5*kn1*dt, t(i-1) + 0.5*dt);
    kn3 = a(xrk(i-1) + 0.5*kn2*dt, t(i-1) + 0.5*dt);
    kn4 = a(xrk(i-1) + kn3*dt, t(i));
    xrk(i) = xrk(i-1) + 1./6*(kn1 + 2*kn2 + 2*kn3 + kn4)*dt;
end

toc

plot(t,x,'linewidth',2,'-r');hold on;
plot(t,xrk,'linewidth',2,'-k');hold off
legend('Euler','4th order RK');
xlabel('time');ylabel('x(t)');

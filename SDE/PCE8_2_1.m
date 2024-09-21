clear all

a = @(x,t) t*x*(2 - x);
dat = @(x,t) x*(2 - x);
dax = @(x,t) t*(2 - 2*x);
ddatt = @(x,t) 0;
ddatx = @(x,t) 2*(1 - x);
ddaxx = @(x,t) -2*t;

tic

n = 16;
t0 = 0.; T = .5;
x0 = 1.;

dt = (T-t0)/n;
x = x0*ones(1,n+1);
x2 = x0*ones(1,n+1);
x3 = x0*ones(1,n+1);
t = t0*ones(1,n+1);

for i = 2:n+1
    xt = x(i-1);
    x(i) = xt + a(xt, t(i-1))*dt;
    xt = x2(i-1);
    x2(i) = xt + a(xt, t(i-1))*dt + ...
           1./2*(dat(xt, t(i-1)) + dax(xt, t(i-1)))*dt*dt;
    xt = x3(i-1);
    x3(i) = xt + a(xt, t(i-1))*dt + ...
           1./2*(dat(xt, t(i-1)) + dax(xt, t(i-1)))*dt*dt + ...
           1./6*(ddatt(xt, t(i-1)) + 2*a(xt, t(i-1))*ddatx(xt, t(i-1)) + ...
           a(xt, t(i-1))*a(xt, t(i-1))*ddaxx(xt, t(i-1)) + ...
           dat(xt, t(i-1))*dax(xt, t(i-1)) + ...
           a(xt, t(i-1))*dax(xt, t(i-1))*dax(xt, t(i-1)))*dt*dt*dt;
    t(i) = t(i-1) + dt;
end

toc

plot(t,x,'linewidth',2,'-r');hold on;
plot(t,x2,'linewidth',2,'-b');
plot(t,x3,'linewidth',2,'-k');hold off
legend('Euler','2nd order','3rd order');
xlabel('time');ylabel('x(t)');

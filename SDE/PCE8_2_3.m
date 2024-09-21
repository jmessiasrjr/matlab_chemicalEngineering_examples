clear all

a = @(x,t)  -3*x + 1;

tic

n = 10;
t0 = 0.; T = 1.;
x0 = 1.;

dt = (T-t0)/n;
x = x0*ones(1,n+1);
xmp = x0*ones(1,n+1);
t = t0*ones(1,n+1);

for i = 2:n+1
    t(i) = t(i-1) + dt;
    x(i) = x(i-1) + a(x(i-1), t(i-1))*dt;
    if (i>2)
        xmp(i) = xmp(i-2) + 2*a(xmp(i-1), t(i-1))*dt;
    else
        xmp(i) = x(i);
    end
end

toc

plot(t,x,'linewidth',2,'-r');hold on;
plot(t,xmp,'linewidth',2,'-k');hold off
legend('Euler','Midpoint');
xlabel('time');ylabel('x(t)');

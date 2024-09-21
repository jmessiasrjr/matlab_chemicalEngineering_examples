clear all

a = @(x,t) -5*x;
Xe = @(x0,t) x0*exp(-5*t);

tic

n = 16;
t0 = 0.; T = 1.;
x0 = 1.;

te = linspace(0,T-t0,1000);
X = Xe(x0,te);
y = X(end);

dt = (T-t0)/n;
x = x0*ones(1,n+1);
xt = x0*ones(1,n+1);
t = t0*ones(1,n+1);

for i = 2:n+1
    x(i) = x(i-1) + a(x(i-1),t(i-1))*dt;
    t(i) = t(i-1) + dt;
    xt(i) = xt(i-1) + 0.5*( a(xt(i-1),t(i-1)) + ...
            a(xt(i-1) + a(xt(i-1),t(i-1))*dt) )*dt;
end

eps = abs(y - x(end));
eps2 = abs(y - xt(end));

fprintf('Global discretization error\n');
fprintf('Euler : %-8.6f #### improved Euler : %-8.6f\n',eps,eps2);

toc

plot(t,x,'linewidth',2,'-r');hold on;
plot(t,xt,'linewidth',2,'-b');
plot(te,X,'linewidth',4,'-k');hold off;
legend('Euler','improved Euler','exact');
xlabel('time');ylabel('x(t)');

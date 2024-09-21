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
xa = x0*ones(1,n+1);
t = t0*ones(1,n+1);

for i = 2:n+1
    x(i) = x(i-1) + a(x(i-1),t(i-1))*dt;
    t(i) = t(i-1) + dt;
    xt(i) = xt(i-1) + 0.5*( a(xt(i-1),t(i-1)) + ...
            a(xt(i-1) + a(xt(i-1),t(i-1))*dt) )*dt;
    if(i>3)
        xa(i) = xa(i-1) + 1./12*( 23*a(xa(i-1),t(i-1)) - ...
                16*a(xa(i-2),t(i-2)) + 5*a(xa(i-3),t(i-3)) )*dt;
    else
        xa(i) = xt(i);
    end
end

eps = abs(y - x(end));
eps2 = abs(y - xt(end));
eps3 = abs(y - xa(end));

fprintf('Global discretization error:\n');
fprintf('Euler       : %-8.6f \n',eps);
fprintf('imp. Euler  : %-8.6f \n',eps2);
fprintf('3-step Adams: %-8.6f \n',eps3);

toc

plot(t,x,'linewidth',2,'-r');hold on;
plot(t,xt,'linewidth',2,'-b');
plot(t,xa,'linewidth',2,'color',[0.4 0.4 0.4],'-p');
plot(te,X,'linewidth',4,'-k');hold off;
legend('Euler','improved Euler','3-step Adams','exact');
xlabel('time');ylabel('x(t)');

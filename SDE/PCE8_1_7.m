clear all

b = 1;

a = @(x,t) -b*x;
Xe = @(x0,t) x0*exp(-b*t);

tic

n = 12;
r = 6; % round off decimal place
t0 = 0.; T = 1.;
x0 = 1.;

te = linspace(0,T-t0,1000);
X = Xe(x0,te);
y = X(end);

q = exp(r*log(10));
x = x0; x2n = x0;
xr = x0; t = t0;

leg = cell(1,n);
leg{1} = 'exact';

fprintf('****************************************\n','#');
fprintf('****ROUND OFF WITH %-d DECIMAL PLACES****\n',r);
fprintf('****************************************\n','#');
figure(1)
plot(te,X,'linewidth',4,'-k');hold on;
fprintf('Global discretization error\n');

for j = 1:n
    m = round(2^j);
    k = 2;
    dt = (T-t0)/m;
    for i = 2:2*m+1
        x2n(i) = round( q*(x2n(i-1) + a(x2n(i-1),t(k-1)/(2-mod(i,2)))*dt/2) )/q;
        if(mod(i,2)==1)
            t(k) = t(k-1) + dt;
            x(k) = round( q*(x(k-1) + a(x(k-1),t(k-1))*dt) )/q;
            xr(k) = round( q*(2*x2n(i) - x(k)) )/q;
            k += 1;
        end
    end
    delta(j) = dt;
    eps(j) = abs(y-x(end));
    eps2(j) = abs(y-xr(end));
    fprintf('n = %-8d : eps = %8.6f\n',m,eps2(j));
    plot(t,x,'linewidth',2,'color',[(n-j)/n 0 j/n]);
    leg{j+1} = num2str(j,'dt = pow(2,-%-d)\n');
    fprintf('Euler scheme   : %-8.6f # ',eps(j));
    fprintf('Romberg interp : %-8.6f \n',eps2(j));
end

hold off;
legend(leg);
xlabel('time');ylabel('x(t)');

toc

figure(2)
plot(log2(delta),log2(eps),'linewidth',2,'-r');hold on;
plot(log2(delta),log2(eps2),'linewidth',2,'-k');hold off;
xlabel('log_2(eps)');ylabel('log_2(dt)');
legend('Euler','Romberg');

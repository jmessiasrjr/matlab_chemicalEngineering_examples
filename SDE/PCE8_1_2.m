clear all

a = @(x,t) -5*x;
Xe = @(x0,t) x0*exp(-5*t);

tic

n = 12;
r = 6; %round off error
t0 = 0.; T = 1.;
x0 = 1.;

te = linspace(0,T-t0,1000);
X = Xe(x0,te);
y = X(end);
q = exp(r*log(10));
t(1) = 0.;x(1) = x0;

leg = cell(1,n);
leg{1} = 'exact';

fprintf('****************************************\n','#');
fprintf('****ROUND OFF WITH %-d DECIMAL PLACES****\n',r);
fprintf('****************************************\n','#');
figure(1)
plot(te,X,'linewidth',4,'-k');hold on;
fprintf('Global discretization error\n');

for i = 1:n
    m = round(2^i);
    for j = 2:m+1
        dt = (T-t0)/(m);
        x(j) = round(q*(x(j-1) + a(x(j-1),t(j-1))*dt))/q;
        t(j) = t(j-1) + dt;
    end
    delta(i) = dt;
    eps(i) = abs(y-x(end));
    fprintf('n = %-8d : eps = %8.6f\n',m,eps(i));
    plot(t,x,'linewidth',2,'color',[(n-i)/n 0 i/n]);
    leg{i+1} = num2str(i,'dt = pow(2,-%-d)\n');
end

hold off;
legend(leg);
xlabel('time');ylabel('x(t)');

toc

figure(2)
plot(log2(delta),log2(eps),'linewidth',2,'-k');
xlabel('log_2(eps)');ylabel('log_2(dt)');

clear all

tic

alpha = 1.5; beta = 1.;
a = @(alpha, x) alpha*x;
b = @(beta, x) beta*x;

X0 = 1.;
t0 = 0.; T = 1;

Xte = @(t, Wt) X0*exp((alpha - 0.5*beta*beta)*t + beta*Wt);

N = 25;
n = 16; % Euler time steps
nex = 256; % exact time steps => must be multiple of n

dt = (T - t0)/n;
t = linspace(t0,T,nex+1);

appstep = round(nex/n);

X = X0*ones(N,nex+1);

DWT = zeros(N,nex);
Xe = X0*ones(N,nex+1);
dte = (T - t0)/nex;

for k = 1:N
    WTe = 0.;
    % Exact solution
    for i = 2:nex+1
        if( mod(i,2)==0 )
            [X1, X2] = GRNnorm;
        else
            X1 = X2;
        end
        DWT(k,i-1) = X1*sqrt(dte);
        WTe += DWT(k,i-1);
        Xe(k,i) = Xte(t(i), WTe);
    end

    % Euler approximation
    for i = 2:n+1
        m = (i-2)*appstep+1;
        X(k,m+appstep) = X(k,m) + a(alpha, X(k,m))*dt + ...
                       b(beta, X(k,m))*sum(DWT(k,m:m+appstep-1));
        if(appstep>1)
            for j = 1:appstep-1 % linear interpolation
                X(k,m+j) = X(k,m) + j/appstep*(X(k,m+appstep) - X(k,m));
            end
        end
    end
    eps(k) = abs(Xe(k,end) - X(k,end));
end

toc

EPS = 1/N*sum(eps);

fprintf('Absolute error : %-10.8f\n',EPS);

plot(t,sum(Xe)/N,'linewidth',2,'-k');hold on;
plot(t,sum(X)/N,'linewidth',2,'-r');hold off;
legend('exact','Euler');
xlabel('time');ylabel('X(t)');


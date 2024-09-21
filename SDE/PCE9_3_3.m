clear all

% May take some time to run ~ 30s average M=20, N=100, nex = 256

tic

alpha = 1.5; beta = .1;
a = @(alpha, x) alpha*x;
b = @(beta, x) beta*x;

X0 = 1.;
t0 = 0.; T = 1;

Xte = @(t, Wt) X0*exp((alpha - 0.5*beta*beta)*t + beta*Wt);

tdist = 1.73; conf = 0.9;
M = 20;
N = 50;
n = 8; % Euler time steps
nex = 64; % exact time steps => must be multiple of n

dt = (T - t0)/n;
t = linspace(t0,T,nex+1);

appstep = round(nex/n);

X = X0*ones(N,nex+1);
Xm = X0*ones(M,nex+1);

DWT = zeros(N,nex);
Xe = X0*ones(N,nex+1);
Xme = X0*ones(M,nex+1);
dte = (T - t0)/nex;

for l = 1:M % batches
    eps = 0;
    for k = 1:N % trajectories
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
        eps += abs(Xe(k,end) - X(k,end));
    end
    Xme(l,:) = sum(Xe)/N;
    Xm(l,:) = sum(X)/N;
    epsj(l) = 1/N*eps;
end

EPS = 1/M*sum(epsj);

var = 0;
for l = 1:M
    var += (epsj(l) - EPS)^2;
end

VAR = 1/(M-1)*var;
Deps = tdist*sqrt(VAR/M);

toc

fprintf('Absolute error : %-10.8f\n',EPS);
fprintf('Variance       : %-10.8f\n',VAR);
fprintf('Interval %-4.1f%% :(%-10.8f %-10.8f)\n',conf*100,EPS-Deps,EPS+Deps)


plot(t,sum(Xme)/M,'linewidth',2,'-k');hold on;
plot(t,sum(Xm)/M,'linewidth',2,'-r');hold off;
legend('exact','Euler');
xlabel('time');ylabel('X(t)');


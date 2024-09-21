clear all

% dXt = -1/2*Xt*dt + Xt*dW1t + Xt*dW2t

% May take some time to run ~ 30s average M=20, N=100, n = 16
% Must be commutative: b_{j1}*db_{j2} = b_{j2}*db_{j1} if j1 != j2

tic

alpha = -.5; beta1 = 1.; beta2 = 1.;
a = @(alpha, x) alpha*x;
b = @(beta, x) beta*x;

X0 = 1.;
t0 = 0.; T = 1;

Xte = @(t, Wt1, Wt2) X0*exp((alpha - 0.5*beta1*beta1 - 0.5*beta2*beta2)*t +...
                     beta1*Wt1 + beta2*Wt2);

tdist = 1.83; conf = 0.9; % see t-student table
M = 10;
N = 50;
n = 32; % Time steps

dt = (T - t0)/n;
t = linspace(t0,T,n+1);

X = X0*ones(N,n+1);
Xm = X0*ones(M,n+1);
Xmil = X0*ones(N,n+1);
Xmmil = X0*ones(M,n+1);

DWT1 = zeros(N,n);
DWT2 = zeros(N,n);
Xe = X0*ones(N,n+1);
Xme = X0*ones(M,n+1);

for l = 1:M % batches
    eps = 0;
    for k = 1:N % trajectories
        WTe1 = 0.; WTe2 = 0.;
        for i = 2:n+1
            [X1, X2] = GRNnorm;
            DWT1(k,i-1) = X1*sqrt(dt);
            DWT2(k,i-1) = X2*sqrt(dt);
            WTe1 += DWT1(k,i-1);
            WTe2 += DWT2(k,i-1);
        % Exact solution
            Xe(k,i) = Xte(t(i), WTe1, WTe2);
        % Euler approximation
            X(k,i) = X(k,i-1) + a(alpha, X(k,i-1))*dt + ...
                           b(beta1, X(k,i-1))*DWT1(k,i-1) + ...
                           b(beta2, X(k,i-1))*DWT2(k,i-1);
        % Milstein scheme
            Xmil(k,i) = Xmil(k,i-1) + a(alpha, Xmil(k,i-1))*dt + ...
                b(beta1, Xmil(k,i-1))*DWT1(k,i-1) + ...
                b(beta2, Xmil(k,i-1))*DWT2(k,i-1) + 0.5*( ...
                b(beta1, Xmil(k,i-1))*beta2*DWT1(k,i-1)*DWT2(k,i-1) + ...
                b(beta2, Xmil(k,i-1))*beta1*DWT2(k,i-1)*DWT1(k,i-1) );
        end
        eps += abs(Xe(k,end) - Xmil(k,end));
    end
    Xme(l,:) = sum(Xe)/N;
    Xm(l,:) = sum(X)/N;
    Xmmil(l,:) = sum(Xmil)/N;
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


plot(t,sum(Xme)/M,'linewidth',4,'-k');hold on;
plot(t,sum(Xm)/M,'linewidth',2,'-r');
plot(t,sum(Xmmil)/M,'linewidth',2,'-b');hold off;
legend('exact','Euler','Milstein');
xlabel('time');ylabel('X(t)');


clear all

% dXt = -1/2*Xt*dt + Xt*dW1t + Xt*dW2t

% May take some time to run ~ 30s average M=20, N=100, n = 16, p = 2

tic

alpha = -.5; beta1 = 1.; beta2 = 0.4;
a = @(alpha, x) alpha*x;
b = @(beta, x) beta*x;

X0 = 1.;
t0 = 0.; T = 1;

Xte = @(t, Wt1, Wt2) X0*exp((alpha - 0.5*beta1*beta1 - 0.5*beta2*beta2)*t +...
                     beta1*Wt1 + beta2*Wt2);
rop = @(sumr2) pi*pi/180 - 1/(2*pi*pi)*sumr2;
Apj = @(r, zeta1, zeta2, ksi1, ksi2, eta1, eta2) 1/r*( zeta1*(sqrt(2)*...
      ksi2 + eta2) - zeta2*(sqrt(2)*ksi1 + eta1) ); 
J12 = @(APJ, ROP, ksi1, ksi2, mu1, mu2, dt) dt/2*(ksi1*ksi2) +...
      dt*sqrt(ROP)*(mu1*ksi2 - mu2*ksi1) + APJ;

tdist = 1.83; conf = 0.9; % see t-student table
M = 10;
N = 50;
n = 16; % Time steps
p = 2; % Truncation index

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
            ROP = 0.; APJ = 0.; sumr2 = 0.;
            [ksi1 ksi2] = GRNnorm;
            [mu1 mu2] = GRNnorm;
            for r = 1:p
                [eta1 eta2] = GRNnorm;
                [zeta1 zeta2] = GRNnorm;
                sumr2 += 1/(r*r);
                APJ += Apj(r, zeta1, zeta2, ksi1, ksi2, eta1, eta2);
            end
            ROP += rop(sumr2);
            J12P = J12(APJ, ROP, ksi1, ksi2, mu1, mu2, dt);
            J21P = DWT1(k,i-1)*DWT2(k,i-1) - J12P;
            I11 = 0.5*(DWT1(k,i-1)*DWT1(k,i-1) - dt);
            I22 = 0.5*(DWT2(k,i-1)*DWT2(k,i-1) - dt);
            Xmil(k,i) = Xmil(k,i-1) + a(alpha, Xmil(k,i-1))*dt + ...
                        b(beta1, Xmil(k,i-1))*DWT1(k,i-1) + ...
                        b(beta2, Xmil(k,i-1))*DWT2(k,i-1) + ...
                        b(beta1, Xmil(k,i-1))*beta1*I11 + ...
                        b(beta1, Xmil(k,i-1))*beta2*J12P + ...
                        b(beta2, Xmil(k,i-1))*beta1*J21P + ...
                        b(beta2, Xmil(k,i-1))*beta2*I22;
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


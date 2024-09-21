clear all

% dXt = a*Xt*dt + b*Xt*dWt

% May take some time to run ~ 30s average M=20, N=100, n = 256

tic

alpha = 1.5; beta = 1.25;
a = @(alpha, x) alpha*x;
b = @(beta, x) beta*x;

X0 = 1.;
t0 = 0.; T = 1;

Xte = @(t, Wt) X0*exp((alpha - 0.5*beta*beta)*t + beta*Wt);

tdist = 1.73; conf = 0.9; % see t-student table
M = 20;
N = 100;
n = 16; % Time steps

dt = (T - t0)/n;
t = linspace(t0,T,n+1);

X = X0*ones(N,n+1);
Xm = X0*ones(M,n+1);
Xmil = X0*ones(N,n+1);
Xmmil = X0*ones(M,n+1);

DWT = zeros(N,n);
Xe = X0*ones(N,n+1);
Xme = X0*ones(M,n+1);

for l = 1:M % batches
    eps = 0;
    for k = 1:N % trajectories
        WTe = 0.;
        for i = 2:n+1
            if(mod(i,2)==0)
                [X1, X2] = GRNnorm;
            else
                X1 = X2;
            end
            DWT(k,i-1) = X1*sqrt(dt);
            WTe += DWT(k,i-1);
        % Exact solution
            Xe(k,i) = Xte(t(i), WTe);
        % Euler Scheme
            X(k,i) = X(k,i-1) + a(alpha, X(k,i-1))*dt + ...
                           b(beta, X(k,i-1))*DWT(k,i-1);
        % Milstein Scheme
            B = b(beta, Xmil(k, i-1)); DB = beta;
            A = a(alpha, Xmil(k,i-1)); DA = alpha;
            [G1, G2] = GRNnorm;
            DZ = sqrt(dt)*dt*(G1+G2/sqrt(3));
            Xmil(k,i) = Xmil(k,i-1) + A*dt + B*DWT(k,i-1) + ...
                           0.5*B*DB*(DWT(k,i-1)*DWT(k,i-1) - dt) + ...
                           0.5*B*DA*DZ + 0.5*(A*DA)*dt*dt + ...
                           A*DB*(DWT(k,i-1)*dt - DZ) + ...
                           0.5*B*(DB*DB)*(1./3*(DWT(k,i-1)*DWT(k,i-1)) -...
                            dt)*DWT(k,i-1);
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
fprintf('Interval %-4.1f%% :(%-10.8f %-10.8f)\n',conf*100,EPS-Deps,EPS+Deps);

figure(1);
plot(t,sum(Xme)/M,'linewidth',4,'-k');hold on;
plot(t,sum(Xm)/M,'linewidth',4,'b');
plot(t,sum(Xmmil)/M,'linewidth',2,'-r');hold off;
title('Averaged values - limit central theorem');
legend('exact','Euler','Milstein');
xlabel('time');ylabel('X(t)');
figure(2);
plot(t,Xe(end,:),'linewidth',2,'-k',t,X(end,:),'linewidth',2,'-b',...
     t,Xmil(end,:),'linewidth',2,'-r');
title('Values of the last simulation');
legend('exact','Euler','Milstein');

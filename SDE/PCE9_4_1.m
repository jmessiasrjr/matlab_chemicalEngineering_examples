clear all

% May take some time to run ~ 30s average M=20, N=100, nex = 256

tic

alpha = 1.5; beta = .1;
a = @(alpha, x) alpha*x;
b = @(beta, x) beta*x;

X0 = 1.;
t0 = 0.; T = 1;

Xte = @(t, Wt) X0*exp((alpha - 0.5*beta*beta)*t + beta*Wt);
meangtx = @(t) X0*exp(alpha*(t-t0));
%gtx = @(x) x;

tdist = 1.68; conf = 0.9; %t-student
M = 40;
N = 50;
n = 32; % Euler time steps

dt = (T - t0)/n;
t = linspace(t0,T,n+1);

X = X0*ones(N,n+1);
Xm = X0*ones(M,n+1);
DWT = zeros(N,n);

for l = 1:M % batches
    mu = 0;
    for k = 1:N % trajectories
        % Euler approximation
        for i = 2:n+1
            if( mod(i,2)==0 )
                [X1, X2] = GRNnorm;
            else
                X1 = X2;
            end
            DWT(k,i-1) = X1*sqrt(dt);
            X(k,i) = X(k,i-1) + a(alpha, X(k,i-1))*dt + ...
                     b(beta, X(k,i-1))*DWT(k,i-1);
        end
        mu += X(k,end) - meangtx(t(end));
    end
    Xm(l,:) = sum(X)/N;
    muj(l) = 1/N*mu;
end

MU = 1/M*sum(muj);

var = 0;
for l = 1:M
    var += (muj(l) - MU)^2;
end

VAR = 1/(M-1)*var;
Dmu = tdist*sqrt(VAR/M);

toc

fprintf('Absolute error : %-10.8f\n',MU);
fprintf('Variance       : %-10.8f\n',VAR);
fprintf('Interval %-4.1f%% :(%-10.8f %-10.8f)\n',conf*100,MU-Dmu,MU+Dmu)

plot(t,sum(Xm)/M,'linewidth',2,'-r');
xlabel('time');ylabel('X(t)');

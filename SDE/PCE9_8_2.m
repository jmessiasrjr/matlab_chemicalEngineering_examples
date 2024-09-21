clear all

tic

alpha = 5.; beta = 1.;
a = @(alpha, x) alpha*x;
b = @(beta, x) beta;

X0 = 1.;
t0 = 0.; T = 1;

Xte = @(t, Wt) X0*exp((alpha)*t) + Wt;

n = 16; % Euler time steps
nex = 256; % exact time steps => must be multiple of n

dt = (T - t0)/n;
t = linspace(t0,T,nex+1);

appstep = round(nex/n);
WTe = 0.;

X = X0*ones(1,nex+1);

DWT = zeros(1,nex);
Xe = X0*ones(1,nex+1);
dte = (T - t0)/nex;
varE = 3./(2.*alpha) + exp(2.*alpha*dte)*(1./(2.*alpha) - dte) + ...
         2.*exp(alpha*dte)*(dte - 1./alpha);
varW = dte;
EEW = dte*(exp(alpha*dte - 1));
Cov = [varE EEW;EEW varW];

% Exact solution
for i = 2:nex+1
    if( mod(i,2)==0 )
        [X1, X2] = GRNcorr(Cov); % Correlated Gaussian random numbers
%        [X1, X2] = GRNnorm; % Indep. Gaussian random numbers 
    else
        X1 = X2;
    end
    DWT(i-1) = X1*sqrt(dte);
    WTe += DWT(i-1);
    Xe(i) = Xte(t(i), WTe);
end

% Euler approximation
for i = 2:n+1
    m = (i-2)*appstep+1;
    X(m+appstep) = X(m) + a(alpha, X(m))*dt + b(beta, X(m))*sum(DWT(m:m+appstep-1));
    if(appstep>1)
        for j = 1:appstep-1 % linear interpolation
            X(m+j) = X(m) + j/appstep*(X(m+appstep) - X(m));
        end
    end
end

toc

plot(t,Xe,'linewidth',2,'-k');hold on;
plot(t,X,'linewidth',2,'-r');hold off;
legend('exact','Euler');
xlabel('time');ylabel('X(t)');

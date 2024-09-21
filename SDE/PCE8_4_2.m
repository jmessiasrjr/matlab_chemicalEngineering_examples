clear all

yt = @(x) x;

tic

n = 4;
L = 1e3;
sint = 40;
s = 4;
xi = 0.4; xf = 0.6;
ti = 0.; tf = 1.;

interv = linspace(-5e-4,5e-4,sint);

q = exp(s*log(10));
dt = (tf - ti)/n;

x0 = linspace(xi, xf, L);
t = linspace(ti, tf, n+1);

y = zeros(1,n+1);
r = zeros(1,L);
N = zeros(1,sint);

for k = 1:L
    y(1) = x0(k);
    for i = 2:n+1
        y(i) = y(i-1) + yt(y(i-1))*dt;
        r(k) += y(i) - round( q*y(i) )/q;
    end
    for j = 1:sint-1
        if( (r(k) > interv(j)) && (r(k) <= interv(j+1)) )
            N(j) += 1;
        end
    end
end

toc

mean = 1/L*sum(r);
xvar = 0.;

for i = 1:L
    xvar += (r(i) - mean)^2;
end

var = 1/(L-1)*xvar;;

fprintf('sample mean : %-10.4e\n',mean);
fprintf('variance    : %-10.4e\n',var);

bar(interv,N);
xlabel('roundoff');ylabel('frequency');

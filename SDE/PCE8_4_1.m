clear all

yk = @(y) pi/3.14*y;

tic

n = 10000;
sint = 40;
s = 4;

interv = linspace(-5*10^(-s-1),5*10^(-s-1),sint);

q = exp(s*log(10));
y0 = 0.1;

y = y0*ones(1,n+1);
r = zeros(1,n+1);
N = zeros(1,sint);

for i = 2:n+1
    y(i) = yk(y(i-1));
    r(i) = y(i) - round( q*y(i) )/q;
    for j = 1:sint-1
        if( (r(i) > interv(j)) && (r(i) <= interv(j+1)) )
            N(j) += 1;
        end
    end
end

R = sum(r);

toc

fprintf('Accumulated roundoff error: %-12.8f\n',R);
figure(1);
bar(interv,N);
xlabel('roundoff');ylabel('frequency');
%figure(2);
%plot(r,'linewidth',2,'-k');
%xlabel('iterations');ylabel('roundoff error');

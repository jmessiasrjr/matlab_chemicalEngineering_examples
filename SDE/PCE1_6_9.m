clear all

P = @(t) [(1+exp(-t))/2 (1-exp(-t))/2; (1-exp(-t))/2 (1+exp(-t))/2];

MA = @(l1,l2) [-l1 l1;l2 -l2];

tic

p = .9;
n = 100;
T = 10;

freq = 0;

x = [1 -2];
po = [p 1-p];
L = [0.5 0.5];%lambda

A = MA(L(1),L(2)); 

for i = 1:n
    X(i) = GRNtp(po(1),x(1),x(2));
    ti = 0.;
    while(ti<T)
        if(X(i)==1)
            dt = GRNexp(A(1,2));
        else
            dt = GRNexp(A(2,1));
        end
        ti += dt;
        if(T-ti>=0) 
            if(X(i)==x(1)) X(i)=x(2); else X(i)=x(1); end 
        end
    end
    if(X(i)==x(1)) freq += 1; end
end
freq /= n;

toc

ya = abs(x(2)-x(1))/10;
plot((1:n),X,'linewidth',2,'-k');
ylim([min(x)-ya, max(x)+ya]);
xlabel('frequency');ylabel('value');

fprintf('Relative frequency = %-10.6f\n',freq);

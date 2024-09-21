function [X,iter] = DILU(A,b,tol,max_iter)

n = length(A);

Ds = zeros(n);
L = zeros(n);
U = zeros(n);

X = zeros(n,1);

res = ones(1,n);

iter = 0;

file=fopen('DILU_iter.dat','w');

for i = 1:n
    Ds(i,i) = A(i,i); %diagonal
end

for i = 1:n
    if(i>1)
        for j = 1:i-1
            L(i,j) = A(i,j); %lower
        end
    end
    if(i<n)
        for j = n:-1:i+1
            U(i,j) = A(i,j); %upper
            if(A(i,j)!=0 && A(j,i)!=0)
                Ds(j,j) -= A(j,i)/A(i,i)*A(i,j);
            end
        end
    end
end

P = (Ds + L)*inv(Ds)*(Ds + U);

rho = max(abs(eig(eye(n) - inv(P)*A))); %spectral radius

if((1-rho)<(1./max_iter))
    fprintf("Spectral radius is high than or equal 1: rho = %10.6f\n",rho);
    break;
else
    fprintf("Spectral radius is: rho = %10.4f\n",rho);
end

while(max(res)>tol)
    iter += 1;
    Xtemp = X;
    fprintf(file,'%4i',iter);
    X = (eye(n) - inv(P)*A)*Xtemp + inv(P)*b;
    for i = 1:n
        res(i) = abs(X(i) - Xtemp(i));
        fprintf(file,'%12.6f',X(i));
    end
    if(iter > max_iter) 
        break;
        fprintf('Stopped after %i iterations',max_iter);
    end
end

fclose(file);

fprintf('Number of iterations: %i\n',iter);

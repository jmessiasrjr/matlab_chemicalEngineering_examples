function [X,iter] = SOR(A,b,omega,tol,max_iter)

n = length(A);

D = zeros(n);
LDomega = zeros(n);
U = zeros(n);
X = zeros(n,1);
res = ones(1,n);

iter = 0;

file=fopen('sor_iter.dat','w');

for i = 1:n
    D(i,i) = A(i,i); %diagonal
    for j = 1:i %lower + diagonal
        if(j<i)
            LDomega(i,j) = A(i,j);
        else
            LDomega(i,j) = 1./omega*A(i,j);
        end
    end
    if(i<n)
        for j = n:-1:i+1
            U(i,j) = A(i,j); %upper
        end
    end
end

LDomega_1 = inv(LDomega); %inverse matrix

rho = max(abs(eig(LDomega_1*(-U+((1./omega-1)*D))))); %spectral radius

if((1-rho)<(1./max_iter))
    fprintf("Spectral radius is high than or equal 1: rho = %10.6f\n",rho);
    break;
else
    fprintf("omega = %8.6f\n",omega);
    fprintf("Spectral radius is: rho = %10.4f\n",rho);
end

while(max(res) > tol)
    iter += 1;
    Xtemp = X;
    fprintf(file,'%4i',iter);
    X = LDomega_1*((-U+((1./omega-1)*D))*Xtemp + b ); % SOR
    for i = 1:n
        res(i) = abs(X(i) - Xtemp(i));
        fprintf(file,'%12.6f',X(i));
    end
    fprintf(file,'\n');
    if(iter > max_iter) 
        break; 
        fprintf('Stopped after %i iterations',max_iter);
    end
end

fclose(file);

fprintf('Number of iterations: %i\n',iter);

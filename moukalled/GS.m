function [X,iter] = GS(A,b,tol,max_iter)

n = length(A);

LD = zeros(n);
U = zeros(n);
X = zeros(n,1);
res = ones(1,n);

iter = 0;

file=fopen('gs_iter.dat','w');

for i = 1:n
    for j = 1:i
        LD(i,j) = A(i,j); %lower + diagonal
    end
    if(i<n)
        for j = n:-1:i+1
            U(i,j) = A(i,j); %upper
        end
    end
end

LD_1 = inv(LD); %inverse matrix

rho = max(abs(eig(-LD_1*U))); %spectral radius

if((1-rho)<(1./max_iter))
    fprintf("Spectral radius is high than or equal 1: rho = %10.6f\n",rho);
    break;
else
    fprintf("Spectral radius is: rho = %10.4f\n",rho);
end

while(max(res) > tol)
    iter += 1;
    Xtemp = X;
    fprintf(file,'%4i',iter);
    X = LD_1*( -U*Xtemp + b ); % Gauss-Seidel
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

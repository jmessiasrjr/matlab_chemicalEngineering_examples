function [X,iter] = Jacobi(A,b,tol,max_iter)

n = length(A);

D = zeros(n);
LU = zeros(n);
X = zeros(n,1);
res = ones(1,n);

iter = 0;

file=fopen('jacobi_iter.dat','w');

for i = 1:n
    D(i,i) = A(i,i); %diagonal
    for j = 1:n
        if(i!=j)
            LU(i,j) = A(i,j); %lower + upper
        end
    end
end

D_1 = inv(D); %inverse matrix

rho = max(abs(eig(-D_1*LU))); %spectral radius

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
    X = D_1*( -LU*Xtemp + b ); % Jacobi
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

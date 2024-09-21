clear all

function phi = F(z, K, f)
    n = length(z);
    phi = 0.;
    for i = 1:n
        phi += z(i)*(1. - K(i))/(1 + f*(K(i) - 1.));
    end
endfunction

X =@(z,f,K) z/(1. + f*(K - 1.));
Y =@(z,f,K) z*K/(1. + f*(K - 1.));

% Values
z = [0.1 0.2 0.3 0.4];
K = [4.2 1.75 0.74 0.34];

% Vaporized fraction
phi = fzero(@(f) F(z, K, f),[0. 1,]);

fprintf("-----------------------------------------\n");
fprintf("PHI : %10.4f\n",phi);
fprintf("-----------------------------------------\n");

for i = 1:length(z)
    x(i) = X(z(i), phi, K(i));
    y(i) = Y(z(i), phi, K(i));
    fprintf("x%-2i : %10.4f\t",i,x(i));
    fprintf("y%-2i : %10.4f\n",i,y(i));
end

fprintf("-----------------------------------------\n");
fprintf("Sum : %10.4f\t",sum(x));
fprintf("      %10.4f\n",sum(y));

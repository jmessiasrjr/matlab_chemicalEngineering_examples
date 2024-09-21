clear all

tic

Q = 100./60;
n = 1.5;
K = 37.7;
k = 0.095;
rhos = 2300;
rhof = 4.43;
mu = 3.5e-5;
Et = 0.85;
u = 15;
I = 1.;
j = 1;

% Gauss-Legendre quadrature
x = [-0.973906 -0.865063 -0.679409 -0.433395 -0.148874 0.148874 0.433395 0.679409 0.865063 0.973906];
w = [0.066671 0.149451 0.219086 0.269267 0.295524 0.295524 0.269267 0.219086 0.149451 0.066671];

X = @(D) 1. - exp(-((D/K)^n));
%DD = @(X) -K*((log(1.-X))^(1./n));
G = @(D,Ds) ((D/Ds)^2.)/(1 + ((D/Ds)^2.));
Dstar = 1.e-7;

while(abs(I(j)-Et)>2.e-3)
    for i=1:length(x)
        h(i) = 1.e-6*fzero(@(D) X(D)-0.5*x(i)-0.5,[0 1000]);
        %h(i) = 1.e-6*real(DD(0.5*x(i)+0.5));
        f(i) = G(h(i),Dstar(j));
        F(i) = w(i)*f(i);
    end
    I(j+1) = 0.5*sum(F);
    Dstar(j+1) = Dstar(j)+1.e-7;
    if(Dstar(j+1)>1.e-4) 
        break; 
    end
    j+=1;
end

P = @(Ds,Dc) k*sqrt(mu/(0.125*u*Dc*(rhos-rhof)))-Ds/Dc;
P2 = @(Ds,Dc,Q) k*sqrt(mu*Dc/(Q*(rhos-rhof)))-Ds/Dc;

Dcyl = fzero(@(Dc) P(Dstar(j),Dc),0.01);
Q1 = u*Dcyl*Dcyl/8;

N = Q/Q1;
Q2 = Q/ceil(N);
Dcyln = fzero(@(Dc) P2(Dstar(j),Dc,Q2),0.01);
unew = 8*Q2/(Dcyln*Dcyln);

toc

fprintf('Number of cyclones\t\t: %-d\n',ceil(N));
fprintf('Cut diameter\t\t\t: %-10.3f µm\n',Dstar(j)*1e6)
fprintf('Diameter of cylindrical section\t: %-10.4f m\n',Dcyl);
fprintf('Velocity\t\t\t: %-10.4f m/s\n',unew);

plot(Dstar,I,'linewidth',2,'-k');
xlabel('Cut diameter');
ylabel('Efficiency');grid;

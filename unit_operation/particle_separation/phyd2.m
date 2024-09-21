clear all

arquivo = fopen('dados.dat','w');

tic

K1 = 21;
n1 = 1.5;
K2 = 3.5;
n2 = 1.2;
rhos1 = 2700;
rhos2 = 2100;
rhof = 995.71;
mu = 7.98e-4;
Dc = 0.05;
DuDc = 0.2; %Du/Dc
k = 0.039;
A = 1.73;
B = 145.;
C = 4.75;
beta = 1200;

DP = 1.01325e5*linspace(1.,4.,4);

Cv1 = 120./rhos1; Cv2 = 45./rhos2;
RL =  B*(DuDc)^C;
Uc = @(delp) sqrt(2*delp/(rhof*beta));

Q = @(uc) pi*uc*Dc*Dc/4.;
X = @(D) 1 - exp(-(D/K)^n);
DS = @(Dc,rhos,Q,f,g) 1e6*Dc*k*sqrt((mu*Dc)/(Q*(rhos-rhof)))*f*g;
fRL = @(rl) 1 + A*rl;
gCV = @(CV) 1/((4.8*(1-CV)^2-3.8*(1-CV))^(0.5));
I = @(Ds,K,n) ((1.13*n/(0.138+n))/(1.44-0.279*n+(K/Ds)))*(K/Ds);

fprintf(arquivo,'  DeltaP\t|  Capacity\t|    D* [µm]\t|   Underflow\t|   Overflow\t| Percentage of M\n');
fprintf(arquivo,'  [Pa]  \t|  flow [m^3/h]\t| M %% \t| A %%\t| M %% \t| A %%\t| M %% \t| A %%\t| of feed in OF %%\n');

for i = 1:length(DP)
    UC = Uc(DP(i));
    Qt(i) = Q(UC);
    frl = fRL(RL);

    gcv = gCV(Cv1+Cv2);
    
    Ds1(i) = DS(Dc,rhos1,Qt(i),frl,gcv);
    Ds2(i) = DS(Dc,rhos2,Qt(i),frl,gcv);

    Et1(i) = I(Ds1(i),K1,n1);
    Et2 = I(Ds2(i),K2,n2);
%    Et1(i) = (1.-RL)*I(Ds1(i),K1,n1)+RL;
%    Et2 = (1.-RL)*I(Ds2(i),K2,n2)+RL;
    Mu(i) = Et1(i)*Cv1/(Et1(i)*Cv1+Et2*Cv2);
    Mo(i) = (1.-Et1(i))*Cv1/((1.-Et1(i))*Cv1+(1.-Et2)*Cv2);
    fprintf(arquivo,'  %-10.4e\t| %-10.4f\t| %-5.2f\t| %-5.2f\t| %-5.2f\t| %-5.2f\t| %-5.2f\t| %-5.2f\t| %-5.2f\n',...
    DP(i),Qt(i)*3.6e3,Ds1(i),Ds2(i),Mu(i)*100,(1.-Mu(i))*100,Mo(i)*100,...
    (1.-Mo(i))*100,(1.-Et1(i))*100);
end

fclose(arquivo);

toc

plot(Ds1,DP,'linewidth',2,'-k',Ds2,DP,'linewidth',2,'-r');
xlabel('Cut diameter [10^{-6} m]');
ylabel('Pressure drop [Pa]');
legend('Dstar M','Dstar A');grid;

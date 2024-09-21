clear all

K = 21.5;
n = 1.35;
rhos1 = 1250;
rhos2 = 2100;
rhof = 1210;
mu = 2.7e-3;
Dc = 0.0508;
DuDc = 0.15; %Du/Dc
kr = 0.039; kb = 0.016;
A = 1.73;
Br = 145.; Bb = 55.3;
Cr = 4.75; Cb = 2.63;
betar = 1200;
betab = 7500;

Cv = 0.05;
DP = 3.1026e5;

V1 = 2./rhos1;V2 = 1./rhos2;
Vt = V1+V2;
x1 = Cv*V1/Vt; x2 = Cv*V2/Vt;

Uc = @(beta) sqrt(2*DP/(rhof*beta));
Q = @(uc) pi*uc*Dc*Dc/4;

X = @(D) 1 - exp(-(D/K)^n);
DS = @(Dc,rhos,Q,k,f,g) 1e6*Dc*k*sqrt((mu*Dc)/(Q*(rhos-rhof)))*f*g;
RL = @(B,C) B*(DuDc)^C;
fRL = @(A,rl) 1 + A*rl;
gCV = @(CV) 1/((4.8*(1-CV)^2-3.8*(1-CV))^(0.5));
I = @(Ds) ((1.13*n/(0.138+n))/(1.44-0.279*n+(K/Ds)))*(K/Ds);

ucr = Uc(betar);
Qr = Q(ucr);
rlr = RL(Br,Cr);
frlr = fRL(A,rlr)

ucb = Uc(betab);
Qb = Q(ucb);
rlb = RL(Bb,Cb);
frlb = fRL(A,rlb)

gcv = gCV(Cv)
gcv1 = gCV(Cv*x1)
gcv2 = gCV(Cv*x2)

Dsr1 = DS(Dc,rhos1,Qr,kr,frlr,gcv)
Dsr2 = DS(Dc,rhos2,Qr,kr,frlr,gcv)
Dsb1 = DS(Dc,rhos1,Qb,kb,frlb,gcv)
Dsb2 = DS(Dc,rhos2,Qb,kb,frlb,gcv)

Etr1 = I(Dsr1)
Etr2 = I(Dsr2)
Etb1 = I(Dsb1)
Etb2 = I(Dsb2)

cur = Etr1*x1/(x1+x2);
aur = Etr2*x2/(x1+x2);

xcur = cur/(cur+aur);
xaur = aur/(cur+aur);

cub = (Etb1)*x1/(x1+x2);
aub = (Etb2)*x2/(x1+x2);

xcub = cub/(cub+aub);
xaub = aub/(cub+aub);

fprintf('Capacity flow (Rietema): \t%-10.4f m^3/h\n',Qr*3.6e3);
fprintf('Capacity flow (Bradley): \t%-10.4f m^3/h\n',Qb*3.6e3);
fprintf('Overflow/feed ratio (Rietema): \t%-10.4f %%\n',(1.-rlr)*100);
fprintf('Overflow/feed ratio (Bradley): \t%-10.4f %%\n',(1.-rlb)*100);
fprintf('Underflow (Rietema) Coal: %-10.4f %%\t Ash: %-10.4f%%\n',xcur*100,xaur*100);
fprintf('Underflow (Bradley) Coal: %-10.4f %%\t Ash: %-10.4f%%\n',xcub*100,xaub*100);

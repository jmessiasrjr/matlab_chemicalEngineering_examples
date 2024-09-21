clear all

g = 981; %[cm]
rhof = 1.15
Dt = 2;

lambda = @(beta,v,Dp) 0.39*exp(6.81*beta)*v/Dp;
Cd = @(v,Dp,rhos) 4./3*(rhos - rhof)*Dp*g/(rhof*v^2);
Re = @(beta,cd) 24*exp(3.54*beta)/((cd^(0.85)-0.43^(0.85))^(1./0.85));
mueff = @(v,Dp,re) Dp*v*rhof/re;
S_l = @(lamb,muef) muef*lamb;

ps = [2.55 2.55 3.98 3.98 7.6 7.6];
D = [.2 .5 .3 .5 .3 .5];
b = [.1 .25 .15 .25 .15 .25];
v = [.72 3.61 3.78 8.87 9.85 22.3];

for i = 1:length(ps)
    l(i) = lambda(b(i),v(i),D(i));
    cd(i) = Cd(v(i),D(i),ps(i));
    re(i) = Re(b(i),cd(i));
    mue(i) = mueff(v(i),D(i),re(i));
    SL(i) = S_l(l(i),mue(i));
    fprintf('lambda = %10.4f | CD     = %10.4f | Re     = %10.4f | mueff  = %10.4f | S(l)   = %10.4f\n',l(i),cd(i),re(i),mue(i),SL(i));
end

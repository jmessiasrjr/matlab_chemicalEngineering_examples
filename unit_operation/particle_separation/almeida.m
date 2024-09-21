clear all

function kp = KP(beta, Reinf)
    if(Reinf < 0.1)
        kp = ((1.-beta)/(1.-0.475*beta))^4;
    elseif(Reinf>1e3)
        kp = 1.-beta^(3./2);
    else
        A = 8.91*exp(2.79*beta);
        B = 1.17e-3 - 0.281*beta;
        kp = 10./(1+A*(Reinf^B));
    end
end

for i = 1:10
    b(i) = i/20;
    for j = 1:1001
        rei(j) = 10^(-1.5+(j-1)/200);
        k(i,j) = KP(b(i), rei(j));
    end
    leg{i} = num2str(b(i),'beta=%-5.2f');
end
semilogx(rei,k,'linewidth',2);
legend(leg,'location','southeast');
xlabel('Re_{inf}');
ylabel('kp'); grid;


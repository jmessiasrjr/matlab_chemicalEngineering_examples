clear all

tic

function uvinf = RZ(eps, Reinf)
    if(Reinf < 0.2)
        uvinf = eps^(3.65);
    elseif(Reinf < 1)
        uvinf = eps^(4.35*Reinf^(-0.03)-1);
    elseif(Reinf < 500)
        uvinf = eps^(4.45*Reinf^(-0.1)-1);
    else
        uvinf = eps^(1.39);
    end
end

function uvinf = PM(eps, Reinf)
    uvinf = eps^(5.93*Reinf^(-0.14));
end

function uvinf = CA(eps, Reinf)
    if(Reinf < 0.2)
        if(eps > 0.9)
            uvinf = 4.8*eps - 3.8;
        else
            uvinf = 0.83*eps^(3.94);
        end
    elseif(Reinf < 500)
        A = 0.28*eps^(-5.96);
        B = 0.35 - 0.33*eps;
        uvinf = 1./(1.+A*Reinf^(-B));
    else
        uvinf = 0.095*exp(2.29*eps);
    end
end

for i = 1:10
    por(i) = 0.5 + (i-1)/20;
    for j = 1:1001
        rei(j) = 10.^(-1.5+(j-1)/200);
        vrz(i,j) = RZ(por(i), rei(j));
        vpm(i,j) = PM(por(i), rei(j));
        vca(i,j) = CA(por(i), rei(j));
    end
    leg{i} = num2str(por(i),'porosity=%-5.2f');
end

toc

semilogx(rei,vrz,'linewidth',2,'-k',rei,vpm,'linewidth',2,'-b',rei,vca,'linewidth',2,'-r');
legend('RZ','PM','CA','location','southeast');
xlabel('Re_{inf}');
ylabel('U/v_{inf}'); grid;

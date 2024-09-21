clear all

function [rop alphap] = smultint(p)

    alphap = 0.; rop = 0.;
    for i = 1:p
        incr = 1/(i*i);
        rop += incr;
        alphap += incr*incr;
    end
    alphap = (pi*pi*pi*pi/90-alphap)/(2*pi*pi);
    rop = pi*pi/6-rop;
    rop = rop/(2*pi*pi);
end

function [J1 J10 J01 J11 J101 J110 J011] = multint(p,rop,alphap,DT,DWT)

    SQDT = sqrt(DT);
    [FI1P MUE1P] = GRNnorm();
    A10 = 0.;B1P = 0.;B11P = 0.;
    for i = 1:p
        [FI1(i) ETA1(i)] = GRNnorm();
        A10 += FI1(i)/i;
        B1P += ETA1(i)/(i*i);
        B11P += (FI1(i)*FI1(i) + ETA1(i)*ETA1(i))/(i*i);
    end
    A10 = -(1/pi)*A10*sqrt(2)*SQDT;
    A10 += -2*SQDT*sqrt(rop)*MUE1P;
    B1P = B1P*SQDT/sqrt(2)+SQDT*sqrt(alphap)*FI1P;
    B11P = B11P/(4*pi*pi);
    J1 = DWT*SQDT;
    J10 = 0.5*DT*(J1+A10);
    J01 = J1*DT-J10;
    J11 = 0.5*J1*J1;
    C11P = 0.;
    R = 0.;
    for i = 1:p
        for j = 1:p
            if(i!=j)
            C11P += (i/(i*i-j*j))*(FI1(i)*FI1(j)/j-ETA1(i)*ETA1(j)*j/i);
            end
        end
    end
    C11P = -C11P/(2.*pi*pi);
    J101 = (DT*DWT)^2/6. - DT*A10*A10/4. + SQDT*DT*DWT*B1P/pi - DT*DT*B11P;
    J110 = (DT*DWT)^2/6. + DT*A10*A10/4. - SQDT*DT*DWT*B1P/(2.*pi) + SQDT*DT*A10*DWT/4. - DT*DT*C11P;
    J011 = J11*DT - J101- J110;
end

n = 100;
T0 = 0.; T = 1;
DT = (T - T0)/n;
p = 10; %truncation index

tic

[rop alphap] = smultint(p);

WT = zeros(1,n+1);
TI = zeros(1,n+1);
J1 = zeros(1,n+1);
J10 = zeros(1,n+1);
J01 = zeros(1,n+1);
J11 = zeros(1,n+1);
J101 = zeros(1,n+1);
J110 = zeros(1,n+1);
J011 = zeros(1,n+1);

for i = 2:n+1
    TI(i) = TI(i-1) + DT;
    if(mod(i,2)==0)
        [DWT1 DWT2] = GRNnorm();
    else
        DWT1 = DWT2;
    end
    [J1(i) J10(i) J01(i) J11(i) J101(i) J110(i) J011(i)] = multint(p,rop,alphap,DT,DWT1);
end

toc

plot(TI,J1,'linewidth',2,'-k');hold on;
plot(TI,J10,'linewidth',2,'-r');
plot(TI,J01,'linewidth',2,'-b');
plot(TI,J11,'linewidth',2,'-g');
plot(TI,J101,'linewidth',2,'.k');
plot(TI,J110,'linewidth',2,'.r');
plot(TI,J011,'linewidth',2,'.b');hold off;
legend('J1','J10','J01','J11','J101','J110','J011');

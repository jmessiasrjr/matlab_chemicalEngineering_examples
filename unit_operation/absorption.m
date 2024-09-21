%  Gas absorber: acetone - air

%          Va,ya
%            ^
%           _|_
%          |   |
% La,xa -->|   |
%          |   |
%          |   |
%          |   |
% Vb,yb -->|   |
%          |___|--> Lb,xb
%
clear all

Vb = 100.; %[mols]
y_abs = 0.97;
r_eq = 1.9; %ratio between line operation and equilibrium data

n = 100; %number of points
g = 3; %polynomial degree approximation

xa = 0.;
xb = 0.1;
yb = 0.3;

Va = (1. - y_abs*yb)*Vb;
Lb = (y_abs*yb*Vb)/xb;
La = Lb - (y_abs*Vb*yb);

ya = (1. - y_abs)/(1./yb - y_abs);

guess = (1. - y_abs)*Vb*yb;

v = linspace(guess, yb*Vb, n);

for i = 1:n
    xnum = v(i) - (1 - y_abs)*Vb*yb;
    x(i) = xnum / (xnum + (1-xb)*Lb);
    y_op(i) = v(i) / (v(i) + (1-yb)*Vb);
    y_eq(i) = y_op(i)/r_eq;
end

q1 = polyfit(x,y_op,g);
q2 = polyfit(x,y_eq,g);

points = linspace(0.,x(end),n);

r1 = polyval(q1, points);
r2 = polyval(q2, points);

plot(points,r1,'-k','linewidth',2,points,r2,'-b','linewidth',2);
xlabel('X');ylabel('Y');grid;
legend('operation line','equilibrium curve');
hold on;

X = 0.; Y = 0.;
i = 2; plates = -1;

while(X < x(end))
    X(i) = X(i-1);
    Y(i) = polyval(q1,X(i));
    X(i+1) = fzero(@(w) Y(i) - polyval(q2,w),0.);
    Y(i+1) = Y(i);
    i += 2;
    plates += 1;
end

plot(X,Y,'-r')
hold off;

fprintf("Number of ideal plates: %6i\n",plates);
fprintf("Entering liquid stream - La: \t%10.4f mols\n",La);
fprintf("Leaving liquid stream - Lb: \t%10.4f mols\n",Lb);
fprintf("Leaving gas stream - Va: \t%10.4f mols\n",Va);

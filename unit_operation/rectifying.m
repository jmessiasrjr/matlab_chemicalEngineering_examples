%  Rectifying benzene - toluene

%           D,yD
%            ^----
%           _|_  |
%          |   |<- 
%          |___|
%  F,zF--> |___|
%          |   |
%          |   |
%          |   |
%          |___|<-
%            |   |
%            v----
%           B,xB 
%
clear all

yeq = @(x, y, alpha) y + (alpha-1)*x*y - alpha*x;

%%% USER DEFINITION %%%%
% Feed stream
F = 200.; %[kgmol/h]
% Molar fraction of benzene
zF = 0.4;
xD = 0.95;
xB = 0.1;

n = 1000; %number of points
g = 8; %polynomial degree approximation to equilibrium curve

% Liquified fraction of feed stream
q = -0.5;
% Relative volatility \alpha_{ij}
alpha = 2.5;
% Reflux ratio
RD = 4.5;

%%%%%%%%%%%%%%%%%%%%%%%%%

% Distillation and bottom streams
D = F*(zF-xB)/(xD-xB);
B = F*(xD-zF)/(xD-xB);

% Line operation of D -- L/D
LOD = [RD/(1.+RD) xD/(1.+RD)];

% Line operation of F
if(abs(q-1.)>1e-4) 
    LOF = [-q/(1.-q) zF/(1.-q)];
else
    LOF = [-1e4 1e4*zF];
end

xi = (LOD(2)-LOF(2))/(LOF(1)-LOD(1));
yi = polyval(LOD,xi);

% Line operation of B -- L/B
LOB = polyfit([xB xi],[xB yi],1);

x = linspace(0.,1.,n);

for i = 1:n
    Yeq(i) = fzero(@(y) yeq(x(i),y,alpha),0.);
end

qeq = polyfit(x,Yeq,g);
qeq2 = qeq;
qeq2(end-1:end) -= LOF;
% Point where LOF and equilibrium curve are crossing by
xmin = fzero(@(w) polyval(qeq2,w),0.);
ymin = polyval(qeq,xmin);
% Calculate minimum L/D
LODmin = polyfit([xmin xD],[ymin xD],1);
RDmin = 1./LODmin(2) - xD;
fprintf("--------------------------------------------------\n");
fprintf("Minimum reflux ratio L/D: \t%10.4f\n",RDmin);
fprintf("--------------------------------------------------\n");

if(RD < RDmin)
    fprintf("Reflux ratio too low!!!:\t%10.4f\n",RD);
    LOD = LODmin;
    LOB = polyfit([xB xmin],[xB ymin],1);
    yi = ymin;
    xi = xmin;
    fprintf("Using minimum Reflux ratio: \t%10.4f \n",RDmin); 
    fprintf("--------------------------------------------------\n");
end

% Graphical solution
plot(x,x,'-k','linewidth',2,x,Yeq,'-b','linewidth',2);
xlabel('X');ylabel('Y');grid;
hold on;
plot([zF xi],[zF yi],'-g','linewidth',2);
plot([xi xD],[yi xD],'-c','linewidth',2);
plot([xB xi],[xB yi],'-m','linewidth',2);

X = xD; Y = xD;
i = 2; plates = 0;

while(X >= xB)
    Y(i) = Y(i-1);
    X(i) = fzero(@(w) Y(i) - polyval(qeq,w),xB);
    if(X(end)<=xi && X(end)>xB)
        Y(i+1) = polyval(LOB,X(i));
    end
    if(X(end)>xi && X(end)>xB)
        Y(i+1) = polyval(LOD,X(i));
    end
    if(X(end)>xB) X(i+1) = X(i); end
    i += 2;
    plates += 1;
    if(plates>1e3) 
        plates = Inf;
        break;
    end
end

plot(X,Y,'-r')
hold off;

fprintf("Number of ideal plates: \t%6i\n",plates);
fprintf("--------------------------------------------------\n");
fprintf("Feed stream - F: \t\t%10.4f kgmol\n",F);
fprintf("Distillation stream - D: \t%10.4f kgmol\n",D);
fprintf("Bottom stream - B: \t\t%10.4f kgmol\n",B);
fprintf("--------------------------------------------------\n");

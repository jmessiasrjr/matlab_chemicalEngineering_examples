% Polar Marsaglia method

function [X1 X2] = GRNnorm()

    W = -1;
    while(W>1 || W<=0)
        V1 = 2*rand()-1;
        V2 = 2*rand()-1;
        W = V1*V1 + V2*V2;
    end
    LW = log(W)/W;
    LW = sqrt(-LW-LW);
    X1 = V1*LW;
    X2 = V2*LW;

end

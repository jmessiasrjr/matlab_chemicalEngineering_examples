function [G1 G2] = GRNcorr(COV)
    A = 0; B = sqrt(COV(1,1));
    D = COV(1,2)/B; C = sqrt(COV(2,2) - D*D);
    [X1 X2] = GRNnorm();
    G1 = A*X1 + B*X2;
    G2 = C*X1 + D*X2;
end

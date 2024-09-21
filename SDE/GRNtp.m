% Generation of two-point random numbers

function [X,rn] = GRNtp(p,x1,x2)
    rn = rand();
    if(rn<=p) 
        X = x1;
    else
        X = x2;
    end
end


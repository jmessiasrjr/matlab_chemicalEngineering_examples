% Generation of exponential random numbers

function [X,rn] = GRNexp(lambda)
    rn = rand();
    if(rn>0) 
        X = -log(rn)/lambda;
    else
        X = 0;
    end
end


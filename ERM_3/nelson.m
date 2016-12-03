function [ r ] = nelson( b, T )
%NELSON interpolate rates and return a yield curve
lambda = 0.7308;
b1 = b(1);
b2 = b(2);
b3 = b(3);

r = b1 + b2 * ((1-exp(-lambda * T))./(lambda * T)) + ...
    b3 * ((1-exp(-lambda * T))./(lambda * T) - exp(-lambda * T));
end



function [ r ] = nelson( b, T )
%NELSON interpolate rates and return a yield curve
lambda = 0.7308;

r = b(1) + b(2) * ((1-exp(-lambda * T))./(lambda * T)) + ...
    b(3) * ((1-exp(-lambda * T))./(lambda * T) - exp(-lambda * T));
end
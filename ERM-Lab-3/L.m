function l = L(x)
    %Alpha, beta, and gamma
    Alpha = 0.000978556;
    Beta = 5.74478E-07;
    Gamma = 0.137829435;
    
    l = exp(-Alpha .* x - Beta .* (exp(Gamma .* x) - 1) ./ Gamma);
end
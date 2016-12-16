function  b  = calibrateNSParameters( T, r)
%CALIBRATENSPARAMETERS This function calibrates the Nelson-Siegel
%parametsers b1, b2, b3.
    % Create matrix of functions f1, f2 and f3
    M = zeros(length(T), 3);
    lambda = 0.7308;
    
    M(:,1) = ones(length(T),1);
    M(:,2) = (1-exp(-lambda * T))./(lambda * T);
    M(:,3) = M(:,2) - exp(-lambda * T)';
    b = ((M' * M) \ M') * r;

end


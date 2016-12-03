function  b  = calibrateNSParameters( T, r)
%CALIBRATENSPARAMETERS This function calibrates the Nelson-Siegel
%parameters b1, b2, b3.

    % Create matrix of functions f1, f2 and f3
    l = length(T);
    M = zeros(l, 3);
    
    M(:,1) = ones(l,1);
    M(:,2) = f2(T');
    M(:,3) = f3(T');
    b = ((M' * M) \ M') * r;


    % Function f2
    function result = f2(T)
        lambda = 0.7308;
        result = (1-exp(-lambda * T))./(lambda * T);
    end
    
    % Function f3
    function result = f3(T)
        lambda = 0.7308;
        result = (1-exp(-lambda * T))./(lambda * T) - exp(-lambda * T);
    end
end



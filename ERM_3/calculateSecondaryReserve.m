function [V, P] = calculateSecondaryReserve(secondaryReserve, portfolioReturn, P, dG, age, my)
        
        tempSecondaryReserve = secondaryReserve * (1 + portfolioReturn);
        riskPremium = my * secondaryReserve;
        
        if age > 65
            payOut = tempSecondaryReserve / (PVUII(age, dG));
            secondaryReserve = tempSecondaryReserve - payOut + riskPremium;
        else
            payOut = 0;
            secondaryReserve = tempSecondaryReserve + P / 12 + riskPremium;
        end
        V = secondaryReserve;
        P = payOut;
end
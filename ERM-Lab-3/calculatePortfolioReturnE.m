function portfolioReturn = calculatePortfolioReturnE(adjScenarioData, rate, BE, firstOrderCashFlow, normalReturns)
%CPPI(1), with bonds matching liability duration in each timestep.
%Prices of assets (Stocks, 5- and 30YR zcb)
prices = adjScenarioData(:,[3:13 end]);

l = 242000;
timeSteps = 121;

%Starting portfolio value
portfolioValue = zeros(l,1);
portfolioValue(1:timeSteps:end,1) = 2200000;

m = 1;

portfolioReturn = zeros(240000,1);

%For i in calculation period
for i = 1:(timeSteps - 1)

    
    %Money to be invested in each stock and the bonds respectively
    stockAll = m .* (portfolioValue(i:timeSteps:end,1) - BE(:,i));
    bondAll = portfolioValue(i:timeSteps:end,1) - stockAll; 
    
    %Timesteps (=monthly 60 years forward)
    t = 1/12:1/12:(60 - (i-1)/12);
    
    %Calculate duration
    rates = squeeze(rate(i,:,1:(end-i+1)));
    durationBE = sum((t .* firstOrderCashFlow((i+1):end)) .* exp(-rates .*t ) ,2)./ BE(:,i);
    
    %Solve equation
    B2 = BE(:,i).*(durationBE-5)./(25*prices(i:121:end, end));
    B1 = (BE(:,i) - B2.*prices(i:121:end ,end))./prices(i:121:end, end-1);

    %Weights for zcb:s (5YR and 30YR)
    w1 = B1 .* prices(i:121:end,end-1)./(B1 .* prices(i:121:end,end-1) + B2 .* prices(i:121:end,end));
    w2 = 1 - w1;

    zcbW = bondAll./(stockAll + bondAll);
    portfolioReturn(i:120:end) = sum([ones(2000,10).*(1-zcbW)/10, zcbW.*w1, zcbW.*w2] .* normalReturns(i:120:end,[1:11 end]),2);
    portfolioValue((i+1):121:end) = portfolioReturn(i:120:end) .* portfolioValue(i:timeSteps:end);
end
portfolioReturn = reshape(portfolioReturn,120, length(portfolioReturn)/120) - 1;
end
function [returns, unAdjMonthlyFinCrisisData] = calculateFinancialCrisisReturn(marketData, strategy)

unAdjFinancialCrisisData = marketData(601:649, :);
unAdjMonthlyFinCrisisData(1, :) = unAdjFinancialCrisisData(1, :);
unAdjMonthlyFinCrisisData(2:13, :) = unAdjFinancialCrisisData(5:4:49, :);

adjMarketData = marketData;
adjMarketData(:, 11:15) = exp(-marketData(:,11:15) .* [5 10 15 20 30]);
financialCrisisData = adjMarketData(601:649, :);
monthlyFinCrisisData = zeros(13, 15);
monthlyFinCrisisData(1, :) = financialCrisisData(1, :);
monthlyFinCrisisData(2:13, :) = financialCrisisData(5:4:49, :);
monthlyFinCrisisLogReturns = diff(log(monthlyFinCrisisData));

finCrisisScenarioPrices = zeros(13,15);
finCrisisScenarioPrices(1, :) = adjMarketData(end, :);
for i = 2:13
    finCrisisScenarioPrices(i, :) = finCrisisScenarioPrices(i-1, :) .* ...
        (1 + monthlyFinCrisisLogReturns(i-1, :));
end
startingPortfolioValue = 2200000;
%Money to be invested in each stock and the bond respectively
if strcmp(strategy, 'BH100')
    stockAll = 0;
    bondAll = startingPortfolioValue/2;
    allocation = [ones(1,10).*stockAll, bondAll]./finCrisisScenarioPrices(1, 13);
    portfolioValue = sum(allocation.*finCrisisScenarioPrices(:,13),2);
    returns = diff(portfolioValue)./portfolioValue(1:end-1,:);
else
    monthlyFinCrisisLogReturns = exp(monthlyFinCrisisLogReturns);

    r = unAdjMonthlyFinCrisisData(:, 11:15);
    zcbReturns = exp(-(r(2:13,:) - r(1:12,:)) .* [5 10 15 20 30]) .* exp(r(2:13,:) * 1/12);
    monthlyFinCrisisLogReturns(:, 11:15) = zcbReturns;
    stockAll = startingPortfolioValue/20;
    bondAll = startingPortfolioValue/2;
    allocation = zeros(13, 11);
    allocation(:,:) = [ones(13,10).*stockAll, ones(13,1)*bondAll];
    
    for T = 2:13
        t = T:13:length(allocation);
        allocation(t,:) = allocation(t-1,:) .* monthlyFinCrisisLogReturns(T-1, 1:11);
    end
    
    portfolioValue = sum(allocation,2);
    returns = diff(portfolioValue)./portfolioValue(1:end-1,:);
    
end
end
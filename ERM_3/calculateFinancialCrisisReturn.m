function returns = calculateFinancialCrisisReturn(marketData, strategy)
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
    stockAll = startingPortfolioValue/20;
    bondAll = startingPortfolioValue/2;
    allocation = [ones(1,10).*stockAll, bondAll]./finCrisisScenarioPrices(1, 1:11);
    portfolioValue = sum(allocation.*finCrisisScenarioPrices(:,1:11),2);
    returns = diff(portfolioValue)./portfolioValue(1:end-1,:);
end
end
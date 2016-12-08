function portfolioReturn = calculatePortfolioReturnA(adjScenarioData)

startingPortfolioValue = 2200000;
%Money to be invested in each stock and the bond respectively
bondAll = startingPortfolioValue;
allocation = bondAll./adjScenarioData(1, 13);
portfolioValue = allocation.*adjScenarioData(:,13);
portfolioValue = reshape(portfolioValue,121,2000);
portfolioReturn = diff(portfolioValue)./portfolioValue(1:end-1,:);
end
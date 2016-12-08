function portfolioReturn = calculatePortfolioReturnB(adjScenarioData)
%Prices of assets
startingPortfolioValue = 2200000;
%Money to be invested in each stock and the bond respectively
stockAll = startingPortfolioValue/20;
bondAll = startingPortfolioValue/2;
allocation = [ones(1,10).*stockAll, bondAll]./adjScenarioData(1, 3:13);
portfolioValue = sum(allocation.*adjScenarioData(:,3:13),2);
portfolioValue = reshape(portfolioValue,121,2000);
portfolioReturn = diff(portfolioValue)./portfolioValue(1:end-1,:);
end

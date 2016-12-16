function portfolioReturn = calculatePortfolioReturnB(normalReturns)
%Prices of assets
startingPortfolioValue = 2200000;
%Money to be invested in each stock and the bond respectively
stockAll = startingPortfolioValue/20;
bondAll = startingPortfolioValue/2;
allocation = zeros(242000, 11);
allocation(1:121:end,:) = [ones(2000,10).*stockAll, ones(2000,1)*bondAll];
for T = 2:121
    t = T:121:length(allocation);
    allocation(t,:) = allocation(t-1,:) .* normalReturns((T-1):120:end,1:11);
end
portfolioValue = sum(allocation,2);
portfolioValue = reshape(portfolioValue,121,2000);
portfolioReturn = diff(portfolioValue)./portfolioValue(1:end-1,:);
end

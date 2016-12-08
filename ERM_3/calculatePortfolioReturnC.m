function portfolioReturn = calculatePortfolioReturnC(adjScenarioData)
%Prices of assets
prices = adjScenarioData(:,3:13);
l = 242000;
timeSteps = 121;

%Starting portfolio value
portfolioValue = zeros(l,1);
portfolioValue(1:timeSteps:end,1) = 2200000;

%Money to be invested in each stock and the bond respectively
stockAll = portfolioValue(1:timeSteps:end,1)./20;
bondAll = portfolioValue(1:timeSteps:end,1)./2;

%Starting allocation
startingPrices = prices(1,:);
allocation = ones(l,11);
allocation(1:timeSteps:end,:) = [ones(1,10).*stockAll, bondAll]./startingPrices;

%For all i in time period
for i = 2:timeSteps
    %Calculate portfolio value
    portfolioValue(i:timeSteps:end) = sum(allocation(i-1:timeSteps:end,:).*prices(i:timeSteps:end,:),2);
    
    %Money to be invested in each stock and the bond respectively
    stockAll = portfolioValue(i:timeSteps:end,1)./20;    
    bondAll =  portfolioValue(i:timeSteps:end,1)./2; 
    
    %Calculate new allocation
    allocation(i:timeSteps:end,:) = horzcat(ones(2000,10).*stockAll, bondAll)./prices(i:timeSteps:end,:);
end

portfolioValue = reshape(portfolioValue,121,2000);
portfolioReturn = diff(portfolioValue)./portfolioValue(1:end-1,:);
end
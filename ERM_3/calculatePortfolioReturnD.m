function portfolioReturn = calculatePortfolioReturnD(adjScenarioData, rate, BE, firstOrderCashFlow)
%CM(50/50) with bond investments matching liability duration in t=0.
disp('Calculating portfolio returns with CM(50/50) duration matching strategy..')
%Prices of assets (Stocks, 5- and 30YR zcb)
prices = adjScenarioData(:,[3:13 end]);

%Timesteps (=monthly 40 years forward)
t = 1/12:1/12:60;
%Calculate duration
rates = zeros(2000,720);
rates(:,:) = rate(1,:,:);
durationBE = sum(ones(2000,720).*(t.*firstOrderCashFlow(2:end)).*exp(-rates.*t),2)./BE(:,1);

%Solve equation
B2 = BE(:,1).*(durationBE-5)./(25*prices(1,end));
B1 = (BE(:,1) - B2.*prices(1,end))./prices(1,end-1);

%Weights for zcb:s (5YR and 30YR)
w1 = B1./(B1+B2);
w2 = 1 - w1;

l = 242000;
timeSteps = 121;

%Starting portfolio value
portfolioValue = zeros(l,1);
portfolioValue(1:timeSteps:end,1) = 2200000;

%Money to be invested in each stock and the bonds respectively
stockAll = portfolioValue(1:timeSteps:end,1)./20;
bondAll = portfolioValue(1:timeSteps:end,1)./2;

%Starting allocation
startingPrices = prices(1,:);
allocation = ones(l,12);
allocation(1:timeSteps:end,:) = [ones(1,10).*stockAll, bondAll.*w1, bondAll.*w2]./startingPrices;

%For i in calculation period
for i = 2:timeSteps
    %Calculate portfolio value
    portfolioValue(i:timeSteps:end) = sum(allocation(i-1:timeSteps:end,:).*prices(i:timeSteps:end,:),2);
    
    
    %Money to be invested in each stock and the bonds respectively
    stockAll = portfolioValue(i:timeSteps:end,1)./20;    
    bondAll =  portfolioValue(i:timeSteps:end,1)./2; 
    
    allocation(i:timeSteps:end,:) = horzcat(ones(2000,10).*stockAll, bondAll.*w1, bondAll.*w2)./prices(i:timeSteps:end,:);
end
%Portfolio value and return (121 rows (timesteps), 2000 columns )
portfolioValue = reshape(portfolioValue,121,2000);
portfolioReturn = diff(portfolioValue)./portfolioValue(1:end-1,:);
end
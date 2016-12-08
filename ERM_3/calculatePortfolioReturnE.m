function portfolioReturn = calculatePortfolioReturnE(adjScenarioData, rate, BE, firstOrderCashFlow)
%CPPI(1), with bonds matching liability duration in each timestep.
%Prices of assets (Stocks, 5- and 30YR zcb)
prices = adjScenarioData(:,[3:13 end]);

%Timesteps (=monthly 40 years forward)
t = 1/12:1/12:60;
%Calculate duration
rates = zeros(2000,720);
rates(:,:) = rate(1,:,:);
%testVar = ones(2000,720).*(t.*firstOrderCashFlow(2:end)).*exp(-rates.*t);
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
m = 1;
stockAll = m .* (portfolioValue(1:timeSteps:end,1)-BE(:,1))./10;
bondAll = portfolioValue(1:timeSteps:end,1) - stockAll.*10;

a = stockAll.*10 + bondAll;

%Starting allocation
startingPrices = prices(1,:);
allocation = ones(l,12);
allocation(1:timeSteps:end,:) = [ones(1,10).*stockAll, bondAll.*w1, bondAll.*w2]./startingPrices;

%For i in calculation period
for i = 2:timeSteps
    %Calculate portfolio value
    portfolioValue(i:timeSteps:end) = sum(allocation(i-1:timeSteps:end,:).*prices(i:timeSteps:end,:),2);
    
    %Money to be invested in each stock and the bonds respectively
    stockAll = m .* (portfolioValue(i:timeSteps:end,1)-BE(:,i))./10;
    bondAll = portfolioValue(i:timeSteps:end,1) - stockAll.*10; 
    
    %Timesteps (=monthly 60 years forward)
    t = 1:(60*12-i+1);
    t = t./12;
    
    %Calculate duration
    rates = zeros(2000,720-(i-1));
    rates(:,:) = rate(i,:,1:(end-i+1));
    durationBE = sum(ones(2000,720-(i-1)).*(t.*firstOrderCashFlow(2:end-(i-1))).*exp(-rates.*t),2)./BE(:,i);
    
    %Solve equation
    B2 = BE(:,i).*(durationBE-5)./(25*prices(i:121:end, end));
    B1 = (BE(:,i) - B2.*prices(i:121:end ,end))./prices(i:121:end, end-1);

    %Weights for zcb:s (5YR and 30YR)
    w1 = B1./(B1+B2);
    w2 = 1 - w1;
    
    allocation(i:timeSteps:end,:) = horzcat(ones(2000,10).*stockAll, bondAll.*w1, bondAll.*w2)./prices(i:timeSteps:end,:);
end
%Portfolio value and return (121 rows (timesteps), 2000 columns )
portfolioValue = reshape(portfolioValue,121,2000);
portfolioReturn = diff(portfolioValue)./portfolioValue(1:end-1,:);
end
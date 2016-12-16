function portfolioReturn = calculatePortfolioReturnD(adjScenarioData, rate, BE, firstOrderCashFlow, normalReturns)
%CM(50/50) with bond investments matching liability duration in t=0.
%Prices of assets (Stocks, 5- and 30YR zcb)
prices = adjScenarioData(:,[3:13 end]);

%Timesteps (=monthly 60 years forward)
t = 1/12:1/12:60;
%Calculate duration
rates(:,:) = squeeze(rate(1,:,:));
durationBE = sum((t.*firstOrderCashFlow(2:end)).*exp(-rates.*t),2)./BE(:,1);

%Solve equation
B2 = BE(:,1).*(durationBE-5)./(25*prices(1,end));
B1 = (BE(:,1) - B2.*prices(1,end))./prices(1,end-1);

%Weights for zcb:s (5YR and 30YR)
w1 = B1*prices(1,end-1)./(B1*prices(1,end-1) + B2*prices(1,end));
w2 = 1 - w1;


%Calculate returns with calculated weights
weights = [ones(1,10)/20 0.5*w1(1) 0.5*w2(1)];
portfolioReturn = sum(normalReturns(:,[1:11 end]) .*  weights,2) - 1;
portfolioReturn = reshape(portfolioReturn,120, length(portfolioReturn)/120);
end
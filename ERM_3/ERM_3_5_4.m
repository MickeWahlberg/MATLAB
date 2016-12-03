%% Comments
% riskpremie + premie -> avkastning -> utbetalning
%% Clear command window
clc
load('data');
load('population')

%If data not already loaded:
if not(exist('popData', 'var'))
    %Load data
    disp('Loading population...')
    f = load('population');
    popData = f.population;
    disp('Population loaded!')
    
end

%If data not already loaded:
if not(exist('marketData', 'var'))
    %Load data
    disp('Loading market data...')
    f = load('MarketData9615');
    marketData = f.MarketData;
    disp('MarketData loaded!')
end

%If data not already loaded:
if not(exist('scenarioData', 'var'))
    %Load data
    disp('Loading scenario data...')
    ScenarioData = load('scenarioData');
    scenarioData = ScenarioData.scenarioData;
    disp('Scenariod data loaded!')
end

%% Parameters
disp('Setting up misc data..')
% Transfer rate and guaranteed return
transferRate = 0.02;
rG = 0.02;

% Spot rate
rS = 0.03;

% Long term asset return
pA = 0.05;

% Starting 1:st and 2:nd order reserves
G0 = 2000000;
V0 = 2200000;

% Log of rates
dG = log(1 + rG);
dS = log(1 + rS);
dV = log(1 + pA);

% Retirement age, payout, and premium
retirementAge = 65;
b = 170619.01695;
P = 120000;

% Makeham parameters
gamma = 0.137829435;
beta = 5.74478E-07;
alfa = 0.000978556;

% Misc data
scenarios = 2000;
periods = 121;
T = [5, 10, 15, 20, 30];

%% First Reserve
t = (60*12:120*12)';

firstReserve = arrayfun(@(x, y) G(x, y), t, ones(length(t),1)*dG);
expectedFirstReserve = firstReserve(:) .* (population(:, 9) + population(:, 10));

%% Transform ScenarioFile
dataFix = scenarioData(:,13:end);
dataFix(dataFix>0.05) = 0.05;
%scenarioData(:,13:end) = dataFix;
adjScenarioData = scenarioData;
adjScenarioData(:,13) = exp(-scenarioData(:,13)*5);
adjScenarioData(:,14) = exp(-scenarioData(:,14)*10);
adjScenarioData(:,15) = exp(-scenarioData(:,15)*15);
adjScenarioData(:,16) = exp(-scenarioData(:,16)*20);
adjScenarioData(:,17) = exp(-scenarioData(:,17)*30);

 

%% Create data for stress test
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
logReturnsFinCrisis = diff(log(finCrisisScenarioPrices));

%% Create returns
a = adjScenarioData(:,3:end);
normalReturns = zeros(240000,15);

for i = 2:121
    t = i:121:242000;
    normalReturns(i-1:120:end,:) = (a(t,:)-a(t-1,:))./a(t-1,:);
end
normalReturns = normalReturns + 1;

%% Assignment B 
%Prices of assets
prices = adjScenarioData(:,3:13);
l = 242000;
timeSteps = 121;

startingPortfolioValue = 2200000;


%Money to be invested in each stock and the bond respectively
stockAll = startingPortfolioValue/20;
bondAll = startingPortfolioValue/2;

allocation = [ones(1,10).*stockAll, bondAll]./adjScenarioData(1, 3:13);

portfolioValue = sum(allocation.*adjScenarioData(:,3:13),2);

portfolioValue = reshape(portfolioValue,121,2000);
portfolioReturnB = diff(portfolioValue)./portfolioValue(1:end-1,:);
%% Assignment C
disp('Calculating portfolio returns with CM(50/50)-strategy..')

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
portfolioReturnC = diff(portfolioValue)./portfolioValue(1:end-1,:);


%% ALM data
secondaryReserve = zeros(scenarios, periods);
tempSecondaryReserve = zeros(scenarios, periods);
riskPremium = zeros(scenarios, periods);
payOut = zeros(scenarios, periods);
expectedSecondaryReserve = zeros(scenarios, periods);
assets = zeros(scenarios, periods);
assetReturns = zeros(scenarios, periods);
vCashFlows = zeros(scenarios, periods);

% Assign start values
secondaryReserve(:, 1) = V0;
tempSecondaryReserve(:, 1) = V0;
expectedSecondaryReserve(:, 1) = V0;
assets(:, 1) = V0;

% Setup empty matrices etc for 5a)
disp('Setting up matrix data..')
rate = zeros(scenarios, 721, 721);
secondaryReserve5a = zeros(scenarios, periods);
tempSecondaryReserve5a = zeros(scenarios, periods);
riskPremium5a = zeros(scenarios, periods);
payOut5a = zeros(scenarios, periods);
expectedSecondaryReserve5a = zeros(scenarios, periods);
assets5a = zeros(scenarios, periods);
assetReturns5a = zeros(scenarios, periods);
vCashFlows5a = zeros(scenarios, periods);
secondaryReserve5a(:, 1) = V0;
tempSecondaryReserve5a(:, 1) = V0;
expectedSecondaryReserve5a(:, 1) = V0;
assets5a(:, 1) = V0;

% Setup matrices etc for 5b)
secondaryReserve5b = zeros(scenarios, periods);
tempSecondaryReserve5b = zeros(scenarios, periods);
riskPremium5b = zeros(scenarios, periods);
payOut5b = zeros(scenarios, periods);
expectedSecondaryReserve5b = zeros(scenarios, periods);
assets5b = zeros(scenarios, periods);
assetReturns5b = zeros(scenarios, periods);
vCashFlows5b = zeros(scenarios, periods);
secondaryReserve5b(:, 1) = V0;
tempSecondaryReserve5b(:, 1) = V0;
expectedSecondaryReserve5b(:, 1) = V0;
assets5b(:, 1) = V0;

% Setup matrices etc for 5c)
secondaryReserve5c = zeros(scenarios, periods);
tempSecondaryReserve5c = zeros(scenarios, periods);
riskPremium5c = zeros(scenarios, periods);
payOut5c = zeros(scenarios, periods);
expectedSecondaryReserve5c = zeros(scenarios, periods);
assets5c = zeros(scenarios, periods);
assetReturns5c = zeros(scenarios, periods);
vCashFlows5c = zeros(scenarios, periods);
secondaryReserve5c(:, 1) = V0;
tempSecondaryReserve5c(:, 1) = V0;
expectedSecondaryReserve5c(:, 1) = V0;
assets5c(:, 1) = V0;

% Setup matrices etc for 5d)
secondaryReserve5d = zeros(scenarios, periods);
tempSecondaryReserve5d = zeros(scenarios, periods);
riskPremium5d = zeros(scenarios, periods);
payOut5d = zeros(scenarios, periods);
expectedSecondaryReserve5d = zeros(scenarios, periods);
assets5d = zeros(scenarios, periods);
assetReturns5d = zeros(scenarios, periods);
vCashFlows5d = zeros(scenarios, periods);
secondaryReserve5d(:, 1) = V0;
tempSecondaryReserve5d(:, 1) = V0;
expectedSecondaryReserve5d(:, 1) = V0;
assets5d(:, 1) = V0;

% Setup matrices etc for 5e)
secondaryReserve5e = zeros(scenarios, periods);
tempSecondaryReserve5e = zeros(scenarios, periods);
riskPremium5e = zeros(scenarios, periods);
payOut5e = zeros(scenarios, periods);
expectedSecondaryReserve5e = zeros(scenarios, periods);
assets5e = zeros(scenarios, periods);
assetReturns5e = zeros(scenarios, periods);
vCashFlows5e = zeros(scenarios, periods);
secondaryReserve5e(:, 1) = V0;
tempSecondaryReserve5e(:, 1) = V0;
expectedSecondaryReserve5e(:, 1) = V0;
assets5e(:, 1) = V0;

% Setup matrices etc for 5f)
secondaryReserve5f = zeros(scenarios, periods);
tempSecondaryReserve5f = zeros(scenarios, periods);
riskPremium5f = zeros(scenarios, periods);
payOut5f = zeros(scenarios, periods);
expectedSecondaryReserve5f = zeros(scenarios, periods);
assets5f = zeros(scenarios, periods);
assetReturns5f = zeros(scenarios, periods);
vCashFlows5f = zeros(scenarios, periods);
secondaryReserve5f(:, 1) = V0;
tempSecondaryReserve5f(:, 1) = V0;
expectedSecondaryReserve5f(:, 1) = V0;
assets5f(:, 1) = V0;

% Setup stress test matrices 5(ii)
tempFinCrisisBh100 = zeros(1,12);
tempFinCrisisBh100(1) = V0;
secondaryFinCrisisBh100 = zeros(1,12);
secondaryFinCrisisBh100(1) = V0;
assetsFinCrisisBh100 = zeros(1,12);
assetsFinCrisisBh100(1) = V0;
rateFinCrisisBh100 = zeros(12,13);

tempFinCrisisBh5050 = zeros(1,12);
tempFinCrisisBh5050(1) = V0;
secondaryFinCrisisBh5050 = zeros(1,12);
secondaryFinCrisisBh5050(1) = V0;
assetsFinCrisisBh5050 = zeros(1,12);
assetsFinCrisisBh5050(1) = V0;
rateFinCrisisBh5050 = zeros(12,13);

bOF5a = zeros(scenarios, periods);
bOF5b = zeros(scenarios, periods);
bOF5c = zeros(scenarios, periods);
bOF5d = zeros(scenarios, periods);
bOF5e = zeros(scenarios, periods);
bOF5f = zeros(scenarios, periods);

%% Nelson-Siegel and Best Estimate
disp('Creating yield curve..')
% Get a 60-year yield curve with nelson-siegel
rate = zeros(121,2000, 720);
for t =  1:121
    for i = 1:2000
        beta = calibrateNSParameters([0.0001 T], ...
            [0.00001 scenarioData(t+((i-1)*121) , 13:17)]');
        rate(t, i, 1:end-(t-1)) = nelson(beta, 1/12 : 1/12 : ((60*12-t+1)/12));
    end
end

% First order cash flow
firstOrderCashFlow = zeros(1, 721);
firstOrderCashFlow(1, 2:end) = b .* population(2:end, 10) ./ 12 + ...
             firstReserve(2:end) .* population(2:end, 8) - ...
             P .* population(2:end, 9) ./ 12;

% Calculate BE
BE = zeros(2000, 121);
for i = 1:121
   t = 1:(720+(1-i));
   rates = zeros(2000,720+(1-i));
   rates(:,:) = rate(i,:,t);
   BE(:,i) = sum(firstOrderCashFlow(t+1).*exp(-rates.*(t./12)),2);
end

%% Assignment D:
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
portfolioReturnD = diff(portfolioValue)./portfolioValue(1:end-1,:);


%% Assignment E:
disp('Calculating portfolio returns with CPPI-strategy..')
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
portfolioReturnE = diff(portfolioValue)./portfolioValue(1:end-1,:);

%% Assignment F:
%The Dart Strategy
%The total capital is invested in a random asset in each time step.
portfolioReturnF = zeros(120, 2000);
asset = randi([1 15], 120,1);
for t = 1:120
    portfolioReturnF(t,:) = normalReturns(t:120:end, asset(t));
end

portfolioReturnF = portfolioReturnF - 1;

%% ERM LABB 4
% Returns ERM_4
    ERM4m0 = 1.0025;
    V0NSS = 1;
    ERM4stocks = marketData(:,1:10);
    ERM4y = ERM4stocks(2:end, :)./ERM4stocks(1:(end-1), :);

    ERM4m = mean(ERM4y(end-100:end,:));

    % Portfolio covariance, variance, and volatility
    ERM4portfolioCov = cov(ERM4y(t-100:t-1,:));
    ERM4portfolioVariance = diag(ERM4portfolioCov);
    ERM4portfolioVolatility = sqrt(ERM4portfolioVariance);

    %Estimate parameters
    ERM4s = length(ERM4portfolioCov);
    ERM4I = ones(ERM4s,1);

    % For same dimensions as in notes
    ERM4m = ERM4m';
    ERM4A = ERM4I' / ERM4portfolioCov * ERM4I;
    ERM4B = ERM4I' / ERM4portfolioCov * ERM4m;
    ERM4C = ERM4m' / ERM4portfolioCov * ERM4m;
    ERM4D = ERM4A * ERM4C - ERM4B^2;
    
    ERM4secondaryReserve = zeros(scenarios, periods);
    ERM4tempSecondaryReserve = zeros(scenarios, periods);
    ERM4riskPremium = zeros(scenarios, periods);
    ERM4payOut = zeros(scenarios, periods);
    ERM4expectedSecondaryReserve = zeros(scenarios, periods);
    ERM4assets = zeros(scenarios, periods);
    ERM4assetReturns = zeros(scenarios, periods);
    ERM4vCashFlows = zeros(scenarios, periods);
    ERM4secondaryReserve(:, 1) = V0;
    ERM4tempSecondaryReserve(:, 1) = V0;
    ERM4expectedSecondaryReserve(:, 1) = V0;
    ERM4assets(:, 1) = V0;
    ERM4bOF = zeros(scenarios, periods);
    ERM4insolvent = zeros(1, scenarios);

    

%% Reserves and assets
disp('Calculating second order reserves and assets...')
insolventA = zeros(1, scenarios);
insolventB = zeros(1, scenarios);
insolventC = zeros(1, scenarios);
insolventD = zeros(1, scenarios);
insolventE = zeros(1, scenarios);
insolventF = zeros(1, scenarios);


for i = 1:scenarios
    for j = 2:periods
        
%         % ERM 4
%         % NO short selling
%         if min(ERM4m) > ERM4m0
%             ERM4weightsNSS(:) = (ERM4m == min(ERM4m).*ERM4secondaryReserve(i, j-1));
%         elseif max(ERM4m) < ERM4m0
%             ERM4weightsNSS(:) = (ERM4m == max(ERM4m).*ERM4secondaryReserve(i, j-1));
%         else
%             ERM4weightsNSS(:) = quadprog(ERM4portfolioCov, [], [], [], [ERM4I';ERM4m'], [ERM4secondaryReserve(i, j-1);ERM4m0*ERM4secondaryReserve(i, j-1)], zeros(ERM4s, 1));
%             ERM4assetsWeightsNSS(:) = quadprog(ERM4portfolioCov, [], [], [], [ERM4I';ERM4m'], [ERM4assets(i, j-1);ERM4m0*ERM4assets(i, j-1)], zeros(ERM4s, 1));
%         end
%         
%         ERM4tempSecondaryReserve(i,j) = sum(ERM4weightsNSS .* ...
%             exp(log(normalReturns((i-1)*120+j-1, 1:10))));
%         ERM4riskPremium(i, j) = data(j, 4) * ERM4secondaryReserve(i, j-1);
%         if j > 61
%             ERM4payOut(i, j) = ERM4tempSecondaryReserve(i, j) / ...
%                 (PVUII(data(j, 2), dG));
%             ERM4secondaryReserve(i, j) = ERM4tempSecondaryReserve(i, j) - ...
%                 ERM4payOut(i, j) + riskPremium(i, j);
%         else
%             ERM4secondaryReserve(i, j) = ERM4tempSecondaryReserve(i, j) - ...
%                 ERM4payOut(i, j) + P / 12 + ERM4riskPremium(i, j);
%         end
%         ERM4expectedSecondaryReserve(i, j) = ERM4secondaryReserve(i, j) * ...
%             (data(j, 9) + data(j, 10));
%         
%         % Calculate value of assets
%         ERM4assetReturns(i, j) = sum(ERM4assetsWeightsNSS .* ...
%             exp(log(normalReturns((i-1)*120+j-1, 1:10)))) - ERM4assets(i, j-1);
%         ERM4vCashFlows(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
%             max(firstReserve(j-1), ERM4secondaryReserve(i, j-1)) - ...
%                 ERM4payOut(i, j) * population(j, 10);
%         ERM4assets(i, j) = ERM4assets(i, j-1) + ERM4vCashFlows(i, j) + ERM4assetReturns(i, j);
%         % Check if insolvent
%         if ERM4assets(i, j) < BE(i, j)
%             ERM4insolvent(i) = 1;
%         end
%         
        
        
        
        
        
        
        % 5(ii)
        if i == 1 && j <= 12
            betaFinCrisisBh100 = calibrateNSParameters([0.001 T], ...
            [0.000001 unAdjMonthlyFinCrisisData((i-1) * 121 + j, 11:15)]');
            rateFinCrisisBh100(j-1, 2:end) = nelson(betaFinCrisisBh100, 1/12 : 1/12 : 1);
        
            % Buy and hold 100% 15 year bond
            tempFinCrisisBh100(j) = secondaryFinCrisisBh100(j-1) * ...
                exp(logReturnsFinCrisis(j, 13));
            riskPremiumFinCrisisBh100(j) = data(j, 4) * secondaryFinCrisisBh100(j-1);
            secondaryFinCrisisBh100(j) = tempFinCrisisBh100(j) - P/12 + ...
                riskPremiumFinCrisisBh100(j);
            
            assetReturnsFinCrisisBh100(i, j) = exp(logReturnsFinCrisis(j, 13)) * ...
            assetsFinCrisisBh100(i, j-1) - assetsFinCrisisBh100(i, j-1);
            vCashFlowsFinCrisisBh100(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
                max(firstReserve(j-1), secondaryFinCrisisBh100(i, j-1));
            assetsFinCrisisBh100(i, j) = assetsFinCrisisBh100(i, j-1) + ...
                vCashFlowsFinCrisisBh100(i, j) + assetReturnsFinCrisisBh100(i, j);
            
            % Buy and hold 50% 5 year bond 50% stocks
            tempFinCrisisBh5050(j) = secondaryFinCrisisBh5050(j-1) * ...
                exp(logReturnsFinCrisis(j, 11) * 0.5 + sum(logReturnsFinCrisis(j, 1:10) * 0.05));
            riskPremiumFinCrisisBh5050(j) = data(j, 4) * secondaryFinCrisisBh5050(j-1);
            secondaryFinCrisisBh5050(j) = tempFinCrisisBh5050(j) - P/12 + ...
                riskPremiumFinCrisisBh5050(j);
            
            assetReturnsFinCrisisBh5050(i, j) = exp(logReturnsFinCrisis(j, 11) * ...
                0.5 + sum(logReturnsFinCrisis(j, 1:10) * 0.05)) * ...
                    assetsFinCrisisBh5050(i, j-1) - assetsFinCrisisBh5050(i, j-1);
            vCashFlowsFinCrisisBh5050(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
                max(firstReserve(j-1), secondaryFinCrisisBh5050(i, j-1));
            assetsFinCrisisBh5050(i, j) = assetsFinCrisisBh5050(i, j-1) + ...
                vCashFlowsFinCrisisBh5050(i, j) + assetReturnsFinCrisisBh5050(i, j);
        end
        
        % 5a) Calculate expected second order reserve (V(t))
        tempSecondaryReserve5a(i,j) = secondaryReserve5a(i, j-1) * ...
            exp(log(normalReturns((i-1)*120+j-1, 13)));
        riskPremium5a(i, j) = data(j, 4) * secondaryReserve5a(i, j-1);
        if j > 61
            payOut5a(i, j) = tempSecondaryReserve5a(i, j) / ...
                (PVUII(data(j, 2), dG));
            secondaryReserve5a(i, j) = tempSecondaryReserve5a(i, j) - ...
                payOut5a(i, j) + riskPremium5a(i, j);
        else
            secondaryReserve5a(i, j) = tempSecondaryReserve5a(i, j) - ...
                payOut5a(i, j) + P / 12 + riskPremium5a(i, j);
        end
        expectedSecondaryReserve5a(i, j) = secondaryReserve5a(i, j) * ...
            (data(j, 9) + data(j, 10));
        
        % 5b) Calculate expected second order reserve (V(t))
        tempSecondaryReserve5b(i,j) = secondaryReserve5b(i, j-1) * ...
            exp(log(1 + portfolioReturnB(j-1, i)));
        riskPremium5b(i, j) = data(j, 4) * secondaryReserve5b(i, j-1);
        
        if j > 61
            payOut5b(i, j) = tempSecondaryReserve5b(i, j) / ...
                (PVUII(data(j, 2), dG));
            secondaryReserve5b(i, j) = tempSecondaryReserve5b(i, j) - ...
                payOut5b(i, j) + riskPremium5b(i, j);
        else
            secondaryReserve5b(i, j) = tempSecondaryReserve5b(i, j) - ...
                payOut5b(i, j) + P / 12 + riskPremium5b(i, j);
        end
        expectedSecondaryReserve5b(i, j) = secondaryReserve5b(i, j) * ...
            (data(j, 9) + data(j, 10));
        
        % 5c) Calculate expected second order reserve (V(t))
        tempSecondaryReserve5c(i,j) = secondaryReserve5c(i, j-1) * ...
            exp(log(1 + portfolioReturnC(j-1, i)));
        riskPremium5c(i, j) = data(j, 4) * secondaryReserve5c(i, j-1);
        
        if j > 61
            payOut5c(i, j) = tempSecondaryReserve5c(i, j) / ...
                (PVUII(data(j, 2), dG));
            secondaryReserve5c(i, j) = tempSecondaryReserve5c(i, j) - ...
                payOut5c(i, j) + riskPremium5c(i, j);
        else
            secondaryReserve5c(i, j) = tempSecondaryReserve5c(i, j) - ...
                payOut5c(i, j) + P / 12 + riskPremium5c(i, j);
        end
        expectedSecondaryReserve5c(i, j) = secondaryReserve5c(i, j) * ...
            (data(j, 9) + data(j, 10));
        
        % 5d) Calculate expected second order reserve (V(t))
        tempSecondaryReserve5d(i,j) = secondaryReserve5d(i, j-1) * ...
            exp(log(1 + portfolioReturnD(j-1, i)));
        riskPremium5d(i, j) = data(j, 4) * secondaryReserve5d(i, j-1);
        
        if j > 61
            payOut5d(i, j) = tempSecondaryReserve5d(i, j) / ...
                (PVUII(data(j, 2), dG));
            secondaryReserve5d(i, j) = tempSecondaryReserve5d(i, j) - ...
                payOut5d(i, j) + riskPremium5d(i, j);
        else
            secondaryReserve5d(i, j) = tempSecondaryReserve5d(i, j) - ...
                payOut5d(i, j) + P / 12 + riskPremium5d(i, j);
        end
        expectedSecondaryReserve5d(i, j) = secondaryReserve5d(i, j) * ...
            (data(j, 9) + data(j, 10));
        
        % 5e) Calculate expected second order reserve (V(t))
        tempSecondaryReserve5e(i,j) = secondaryReserve5e(i, j-1) * ...
            (1 + portfolioReturnE(j-1, i));
        riskPremium5e(i, j) = data(j, 4) * secondaryReserve5e(i, j-1);
        
        if j > 61
            payOut5e(i, j) = tempSecondaryReserve5e(i, j) / ...
                (PVUII(data(j, 2), dG));
            secondaryReserve5e(i, j) = tempSecondaryReserve5e(i, j) - ...
                payOut5e(i, j) + riskPremium5e(i, j);
        else
            secondaryReserve5e(i, j) = tempSecondaryReserve5e(i, j) - ...
                payOut5e(i, j) + P / 12 + riskPremium5e(i, j);
        end
        expectedSecondaryReserve5e(i, j) = secondaryReserve5e(i, j) * ...
            (data(j, 9) + data(j, 10));
        
        % 5f) Calculate expected second order reserve (V(t))
        tempSecondaryReserve5f(i,j) = secondaryReserve5f(i, j-1) * ...
            (1 + portfolioReturnF(j-1, i));
        riskPremium5f(i, j) = data(j, 4) * secondaryReserve5f(i, j-1);
        
        if j > 61
            payOut5f(i, j) = tempSecondaryReserve5f(i, j) / ...
                (PVUII(data(j, 2), dG));
            secondaryReserve5f(i, j) = tempSecondaryReserve5f(i, j) - ...
                payOut5f(i, j) + riskPremium5f(i, j);
        else
            secondaryReserve5f(i, j) = tempSecondaryReserve5f(i, j) - ...
                payOut5f(i, j) + P / 12 + riskPremium5f(i, j);
        end
        expectedSecondaryReserve5f(i, j) = secondaryReserve5f(i, j) * ...
            (data(j, 9) + data(j, 10));
        
        % 5a) Calculate value of assets
        assetReturns5a(i, j) = exp(log(normalReturns((i-1)*120+j-1, 13))) * ...
            assets5a(i, j-1) - assets5a(i, j-1);
        vCashFlows5a(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
            max(firstReserve(j-1), secondaryReserve5a(i, j-1)) - ...
                payOut5a(i, j) * population(j, 10);
        assets5a(i, j) = assets5a(i, j-1) + vCashFlows5a(i, j) + assetReturns5a(i, j);
        % Check if insolvent
        if assets5a(i, j) < BE(i, j)
            insolventA(i) = 1;
        end
        
        % 5b) Calculate value of assets
        assetReturns5b(i, j) = exp(log(1 + portfolioReturnB(j-1, i))) * ...
                    assets5b(i, j-1) - assets5b(i, j-1);
        vCashFlows5b(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
            max(firstReserve(j-1), secondaryReserve5b(i, j-1)) - ...
                payOut5b(i, j) * population(j, 10);
        assets5b(i, j) = assets5b(i, j-1) + vCashFlows5b(i, j) + assetReturns5b(i, j);
        % Check if insolvent
        if assets5b(i, j) < BE(i, j)
            insolventB(i) = 1;
        end
        
        % 5c) Calculate value of assets
        assetReturns5c(i, j) = exp(log(1 + portfolioReturnC(j-1, i))) * ...
                    assets5c(i, j-1) - assets5c(i, j-1);
        vCashFlows5c(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
            max(firstReserve(j-1), secondaryReserve5c(i, j-1)) - ...
                payOut5c(i, j) * population(j, 10);
        assets5c(i, j) = assets5c(i, j-1) + vCashFlows5c(i, j) + assetReturns5c(i, j);
        % Check if insolvent
        if assets5c(i, j) < BE(i, j)
            insolventC(i) = 1;
        end
        
        % 5d) Calculate value of assets
        assetReturns5d(i, j) = exp(log(1 + portfolioReturnD(j-1, i))) * ...
                    assets5d(i, j-1) - assets5d(i, j-1);
        vCashFlows5d(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
            max(firstReserve(j-1), secondaryReserve5d(i, j-1)) - ...
                payOut5d(i, j) * population(j, 10);
        assets5d(i, j) = assets5d(i, j-1) + vCashFlows5d(i, j) + assetReturns5d(i, j);
        % Check if insolvent
        if assets5d(i, j) < BE(i, j)
            insolventD(i) = 1;
        end
        
        % 5e) Calculate value of assets
        assetReturns5e(i, j) = (1 + portfolioReturnE(j-1, i)) * ...
                    assets5e(i, j-1) - assets5e(i, j-1);
        vCashFlows5e(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
            max(firstReserve(j-1), secondaryReserve5e(i, j-1)) - ...
                payOut5e(i, j) * population(j, 10);
        assets5e(i, j) = assets5e(i, j-1) + vCashFlows5e(i, j) + assetReturns5e(i, j);
        % Check if insolvent
        if assets5e(i, j) < BE(i, j)
            insolventE(i) = 1;
        end
        
        % 5f) Calculate value of assets
        assetReturns5f(i, j) = (1 + portfolioReturnF(j-1, i)) * ...
                    assets5f(i, j-1) - assets5f(i, j-1);
        vCashFlows5f(i, j) = (P / 12) * population(j, 9) - population(j, 8) * ...
            max(firstReserve(j-1), secondaryReserve5f(i, j-1)) - ...
                payOut5f(i, j) * population(j, 10);
        assets5f(i, j) = assets5f(i, j-1) + vCashFlows5f(i, j) + assetReturns5f(i, j);
        % Check if insolvent
        if assets5f(i, j) < BE(i, j)
            insolventF(i) = 1;
        end
        
        if j <= 121
            bOF5a(i, j) = assets5a(i, j) - BE(i, j);
            bOF5b(i, j) = assets5b(i, j) - BE(i, j);
            bOF5c(i, j) = assets5b(i, j) - BE(i, j);
            bOF5d(i, j) = assets5c(i, j) - BE(i, j);
            bOF5e(i, j) = assets5d(i, j) - BE(i, j);
            bOF5f(i, j) = assets5e(i, j) - BE(i, j);
            ERM4bOF(i, j) = ERM4assets(i, j) - BE(i, j);
        end
        
    end
end


nrOfTimesInsolventA = sum(insolventA(:) == 1);
nrOfTimesInsolventB = sum(insolventB(:) == 1);
nrOfTimesInsolventC = sum(insolventC(:) == 1);
nrOfTimesInsolventD = sum(insolventD(:) == 1);
nrOfTimesInsolventE = sum(insolventE(:) == 1);
nrOfTimesInsolventF = sum(insolventF(:) == 1);

insolvenceRatioA = nrOfTimesInsolventA / (scenarios);
insolvenceRatioB = nrOfTimesInsolventB / (scenarios);
insolvenceRatioC = nrOfTimesInsolventC / (scenarios);
insolvenceRatioD = nrOfTimesInsolventD / (scenarios);
insolvenceRatioE = nrOfTimesInsolventE / (scenarios);
insolvenceRatioF = nrOfTimesInsolventF / (scenarios);

%% ERM 4
ERM4nrOfTimesInsolvent = sum(ERM4insolvent(:) == 1);
ERM4insolvenceRatio = ERM4nrOfTimesInsolvent / (scenarios);
ERM4bOFYr = zeros(scenarios, 11);
ERM4bOFYr(:, 1) = ERM4assets(:, 1) - BE(:, 1);
ERM4bOFYr(:, 2:end) = ERM4assets(:, 13:12:121) - BE(:, 13:12:121);
ERM4deltaBOFYr = diff(ERM4bOFYr')';
ERM4sCRYr = reshape(ERM4deltaBOFYr, 1, length(ERM4deltaBOFYr(1,:)) * length(ERM4deltaBOFYr(:,1)));
ERM4sCRYr = sort(ERM4sCRYr, 'descend')';
ERM4valueAtRiskSCRYr = ERM4sCRYr(length(ERM4sCRYr) * 0.95);

%% BoF calculation
disp('Calculating BoF and SCR...')
% Calculate BoF per year (not sure which to use yet)
bOF5aYr = zeros(scenarios, 11);
bOF5aYr(:, 1) = assets5a(:, 1) - BE(:, 1);
bOF5aYr(:, 2:end) = assets5a(:, 13:12:121) - BE(:, 13:12:121);
bOF5bYr = zeros(scenarios, 11);
bOF5bYr(:, 1) = assets5b(:, 1) - BE(:, 1);
bOF5bYr(:, 2:end) = assets5b(:, 13:12:121) - BE(:, 13:12:121);
bOF5cYr = zeros(scenarios, 11);
bOF5cYr(:, 1) = assets5c(:, 1) - BE(:, 1);
bOF5cYr(:, 2:end) = assets5c(:, 13:12:121) - BE(:, 13:12:121);
bOF5dYr = zeros(scenarios, 11);
bOF5dYr(:, 1) = assets5d(:, 1) - BE(:, 1);
bOF5dYr(:, 2:end) = assets5d(:, 13:12:121) - BE(:, 13:12:121);
bOF5eYr = zeros(scenarios, 11);
bOF5eYr(:, 1) = assets5e(:, 1) - BE(:, 1);
bOF5eYr(:, 2:end) = assets5e(:, 13:12:121) - BE(:, 13:12:121);
bOF5fYr = zeros(scenarios, 11);
bOF5fYr(:, 1) = assets5f(:, 1) - BE(:, 1);
bOF5fYr(:, 2:end) = assets5f(:, 13:12:121) - BE(:, 13:12:121);
bOFFinCrisisBh100 = assetsFinCrisisBh100(1:11:12) - BE(1, 1:11:12);
bOFFinCrisisBh5050 = assetsFinCrisisBh5050(1:11:12) - BE(1, 1:11:12);

% Delta BoF and SCR for strategy a)
deltaBOF5aYr = diff(bOF5aYr')';
sCR5aYr = reshape(deltaBOF5aYr, 1, length(deltaBOF5aYr(1,:)) * length(deltaBOF5aYr(:,1)));
sCR5aYr = sort(sCR5aYr, 'descend')';
valueAtRiskSCR5aYr = sCR5aYr(length(sCR5aYr) * 0.95);

% Delta BoF and SCR for strategy b)
deltaBOF5bYr = diff(bOF5bYr')';
sCR5bYr = reshape(deltaBOF5bYr, 1, length(deltaBOF5bYr(1,:)) * length(deltaBOF5bYr(:,1)));
sCR5bYr = sort(sCR5bYr, 'descend')';
valueAtRiskSCR5bYr = sCR5bYr(length(sCR5bYr) * 0.95);

% Delta BoF and SCR for strategy c)
deltaBOF5cYr = diff(bOF5cYr')';
sCR5cYr = reshape(deltaBOF5cYr, 1, length(deltaBOF5cYr(1,:)) * length(deltaBOF5cYr(:,1)));
sCR5cYr = sort(sCR5cYr, 'descend')';
valueAtRiskSCR5cYr = sCR5cYr(length(sCR5cYr) * 0.95);

% Delta BoF and SCR for strategy d)
deltaBOF5dYr = diff(bOF5dYr')';
sCR5dYr = reshape(deltaBOF5dYr, 1, length(deltaBOF5dYr(1,:)) * length(deltaBOF5dYr(:,1)));
sCR5dYr = sort(sCR5dYr, 'descend')';
valueAtRiskSCR5dYr = sCR5dYr(length(sCR5dYr) * 0.95);

% Delta BoF and SCR for strategy e)
deltaBOF5eYr = diff(bOF5eYr')';
sCR5eYr = reshape(deltaBOF5eYr, 1, length(deltaBOF5eYr(1,:)) * length(deltaBOF5eYr(:,1)));
sCR5eYr = sort(sCR5eYr, 'descend')';
valueAtRiskSCR5eYr = sCR5eYr(length(sCR5eYr) * 0.95);

% Delta BoF and SCR for strategy f)
deltaBOF5fYr = diff(bOF5fYr')';
sCR5fYr = reshape(deltaBOF5fYr, 1, length(deltaBOF5fYr(1,:)) * length(deltaBOF5fYr(:,1)));
sCR5fYr = sort(sCR5fYr, 'descend')';
valueAtRiskSCR5fYr = sCR5fYr(length(sCR5fYr) * 0.95);

% Delta BoF during financial crisis
deltaBOFFinCrisisBh100 = diff(bOFFinCrisisBh100);
deltaBOFFinCrisisBh5050 = diff(bOFFinCrisisBh5050);

%% Mean stuff
bOF5aYr(find(((bOF5aYr==inf)+(bOF5aYr==-inf)+isnan(bOF5aYr))==1)) = 0;
bOF5bYr(find(((bOF5bYr==inf)+(bOF5bYr==-inf)+isnan(bOF5bYr))==1)) = 0;
bOF5cYr(find(((bOF5cYr==inf)+(bOF5cYr==-inf)+isnan(bOF5cYr))==1)) = 0;
bOF5dYr(find(((bOF5dYr==inf)+(bOF5dYr==-inf)+isnan(bOF5dYr))==1)) = 0;
bOF5eYr(find(((bOF5eYr==inf)+(bOF5eYr==-inf)+isnan(bOF5eYr))==1)) = 0;
bOF5fYr(find(((bOF5fYr==inf)+(bOF5fYr==-inf)+isnan(bOF5fYr))==1)) = 0;



meanA = mean(mean(insolventA .* bOF5aYr(:,end)));
meanB = mean(mean(insolventB .* bOF5bYr(:,end)));
meanC = mean(mean(insolventC .* bOF5cYr(:,end)));
meanD = mean(mean(insolventD .* bOF5dYr(:,end)));
meanE = mean(mean(insolventE .* bOF5eYr(:,end)));
meanF = mean(mean(insolventF .* bOF5fYr(:,end)));

%% Plot
figure
scatter(insolvenceRatioA, meanA)
hold on
scatter(insolvenceRatioB, meanB)
scatter(insolvenceRatioC, meanC)
scatter(insolvenceRatioD, meanD)
scatter(insolvenceRatioE, meanE)
scatter(insolvenceRatioF, meanF)
title('Mean-insolvency probability diagram');
ylabel('Mean of own funds');
xlabel('Insolvency probability');
legend('BH(100/0','BH(50/50)', 'CM(50/50)', 'CH(50/50) with duration matching', 'CPPI(1)', 'The Dart throwing strategy');

%% Plot more
figure
plot(expectedSecondaryReserve5a(1,:), '-');
hold on
plot(assets5a(1,:));
plot(BE(1, 1:121));
plot(bOF5a(1, 2:end));
plot(firstReserve(1:121));
title('Strategy BH(100/0) - scenario 1');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5a(2,:), '-');
hold on
plot(assets5a(2,:));
plot(BE(2, 1:121));
plot(bOF5a(2, 2:end));
plot(firstReserve(1:121));
title('Strategy BH(100/0) - scenario 2');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5b(1,:), '-');
hold on
plot(assets5b(1,:));
plot(BE(1, 1:121));
plot(bOF5b(1, 2:end));
plot(firstReserve(1:121));
title('Strategy BH(50/50) - scenario 1');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5b(2,:), '-');
hold on
plot(assets5b(2,:));
plot(BE(2, 1:121));
plot(bOF5b(2, 2:end));
plot(firstReserve(1:121));
title('Strategy BH(50/50) - scenario 2');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5c(1,:), '-');
hold on
plot(assets5c(1,:));
plot(BE(1, 1:121));
plot(bOF5c(1, 2:end));
plot(firstReserve(1:121));
title('Strategy CM(50/50) - scenario 1');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5c(2,:), '-');
hold on
plot(assets5c(2,:));
plot(BE(2, 1:121));
plot(bOF5c(2, 2:end));
plot(firstReserve(1:121));
title('Strategy CM(50/50) - scenario 2');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5d(1,:), '-');
hold on
plot(assets5d(1,:));
plot(BE(1, 1:121));
plot(bOF5d(1, 2:end));
plot(firstReserve(1:121));
title('Strategy CM(50/50) with duration matching - scenario 1');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5d(2,:), '-');
hold on
plot(assets5d(2,:));
plot(BE(2, 1:121));
plot(bOF5d(2, 2:end));
plot(firstReserve(1:121));
title('Strategy CM(50/50) with duration matching - scenario 2');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');


figure
plot(expectedSecondaryReserve5e(1,:), '-');
hold on
plot(assets5e(1,:));
plot(BE(1, 1:121));
plot(bOF5e(1, 2:end));
plot(firstReserve(1:121));
title('Strategy CPPI(1) with duration matching - scenario 1');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5e(2,:), '-');
hold on
plot(assets5e(2,:));
plot(BE(2, 1:121));
plot(bOF5e(2, 2:end));
plot(firstReserve(1:121));
title('Strategy CPPI(1) with duration matching - scenario 2');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5f(1,:), '-');
hold on
plot(assets5f(1,:));
plot(BE(1, 1:121));
plot(bOF5f(1, 2:end));
plot(firstReserve(1:121));
title('Dart Throwing Strategy - scenario 1');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

figure
plot(expectedSecondaryReserve5f(2,:), '-');
hold on
plot(assets5f(2,:));
plot(BE(2, 1:121));
plot(bOF5f(2, 2:end));
plot(firstReserve(1:121));
title('Dart Throwing Strategy - scenario 2');
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');

disp('Finished!')



%% Functions

function l = L(x)
    %Alpha, beta, and gamma
    Alpha = 0.000978556;
    Beta = 5.74478E-07;
    Gamma = 0.137829435;
    
    l = exp(-Alpha .* x - Beta .* (exp(Gamma .* x) - 1) ./ Gamma);
end

%Returns present value of premiums (1 kr)
function result = PVP(x, dG)
    %Retirement and time to retirement
    z = 65*12;
    u = z - x;
    
        t = 1: u;
        result = (u>0).*sum((L((x + t)./12) ./ L(x/12)) .* exp(-dG .* t/12));
end

function result = PVU(x, dG)
    %Retirement and time to retirement
    z = 65*12;
    u = z - x;
    S = 12*(120 - 65);

    t = (u*(u>0) + 1):(u + S);
    result = sum((L((x + t)./12) ./ L(x/12)) .* exp(-dG .* t/12));
end

function result = PVUII(x, dG)
    %Retirement and time to retirement
    z = 65;
    u = z - x;
    S = 120 - 65;
   
    
    if x > 65
        t = 0:(1/12): u + S;
        result = sum((L(x + t) ./ L(x)) .* exp(-dG * t));
    else 
        result = 0;
    end
      
end

%First order reserve for year x
function firstOrderReserve = G(x, dG)
    b = 170619.01695;
    P = 120000;
    firstOrderReserve = (b / 12) .* PVU(x, dG) - (P / 12) .* PVP(x,dG);
end

function result = q(x)
    result = 1 - (L(x+1/12)/L(x));
end



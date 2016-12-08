%% Comments
% riskpremie + premie -> avkastning -> utbetalning
% ratios: 5%, 25%, 24% 1% ??
% scr: 38k, 250k, 240k, ....
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
transferRate = 0.02;
rG = 0.02;
rS = 0.03;
pA = 0.05;
G0 = 2000000;
V0 = 2200000;
dG = log(1 + rG);
dS = log(1 + rS);
dV = log(1 + pA);
retirementAge = 65;
b = 170619.01695;
P = 120000;
gamma = 0.137829435;
beta = 5.74478E-07;
alfa = 0.000978556;
scenarios = 2000;
periods = 121;
T = [5, 10, 15, 20, 30];

%% First Reserve
t = (60*12:120*12)';
firstReserve = arrayfun(@(x, y) G(x, y), t, ones(length(t),1)*dG);
expectedFirstReserve = firstReserve(:) .* (population(:, 9) + population(:, 10));

%% Transform ScenarioFile
adjScenarioData = scenarioData;
adjScenarioData(:,13) = exp(-scenarioData(:,13)*5);
adjScenarioData(:,14) = exp(-scenarioData(:,14)*10);
adjScenarioData(:,15) = exp(-scenarioData(:,15)*15);
adjScenarioData(:,16) = exp(-scenarioData(:,16)*20);
adjScenarioData(:,17) = exp(-scenarioData(:,17)*30);

%% Create returns
a = adjScenarioData(:,3:end);
normalReturns = zeros(240000,15);

for i = 2:121
    t = i:121:242000;
    normalReturns(i-1:120:end,:) = (a(t,:)-a(t-1,:))./a(t-1,:);
end
normalReturns = normalReturns + 1;

%% Nelson-Siegel and Best Estimate
disp('Creating yield curve..')
% Get a 60-year yield curve with nelson-siegel
rate = zeros(121,2000, 720);
for t =  1:121
    for i = 1:2000
        beta = calibrateNSParameters([T], ...
            [scenarioData(t+((i-1)*121) , 13:17)]');
        rate(t, i, 1:end-(t-1)) = nelson(beta, 1/12 : 1/12 : ((60*12-t+1)/12));
    end
end

% First order cash flow
firstOrderCashFlow = zeros(1, 721);
firstOrderCashFlow(1, 2:end) = b .* population(2:end, 10) ./ 12 + firstReserve(2:end) .* ...
    population(2:end, 8) - P .* population(2:end, 9) ./ 12;

% Calculate BE
BE = zeros(2000, 121);
for i = 1:121
   t = 1:(720+(1-i));
   rates = zeros(2000,720+(1-i));
   rates(:,:) = rate(i,:,t);
   BE(:,i) = sum(firstOrderCashFlow(t+1).*exp(-rates.*(t./12)),2);
end

%% Assignment A
portfolioReturnA = calculatePortfolioReturnA(adjScenarioData);

%% Assignment B 
portfolioReturnB = calculatePortfolioReturnB(adjScenarioData);

%% Assignment C
disp('Calculating portfolio returns with CM(50/50)-strategy..')
portfolioReturnC = calculatePortfolioReturnC(adjScenarioData);

%% Assignment D:
%CM(50/50) with bond investments matching liability duration in t=0.
disp('Calculating portfolio returns with CM(50/50) duration matching strategy..')
portfolioReturnD = calculatePortfolioReturnD(adjScenarioData, rate, BE, firstOrderCashFlow);

%% Assignment E:
portfolioReturnE = calculatePortfolioReturnE(adjScenarioData, rate, BE, firstOrderCashFlow);

%% Assignment F:
%The Dart Strategy
portfolioReturnF = calculatePortfolioReturnF(normalReturns);

%% Create returns for stress test
portfolioReturnFinCrisisA = calculateFinancialCrisisReturn(marketData, 'BH100');
portfolioReturnFinCrisisB = calculateFinancialCrisisReturn(marketData, 'BH5050');

%% ALM data
nStrategies = 6;
secondaryReserve = zeros(nStrategies, scenarios, periods);
tempSecondaryReserve = zeros(nStrategies, scenarios, periods);
riskPremium = zeros(nStrategies, scenarios, periods);
payOut = zeros(nStrategies, scenarios, periods);
expectedSecondaryReserve = zeros(nStrategies, scenarios, periods);
assets = zeros(nStrategies, scenarios, periods);
assetReturns = zeros(nStrategies, scenarios, periods);
vCashFlows = zeros(nStrategies, scenarios, periods);
bOF = zeros(nStrategies, scenarios, periods);

% Assign start values
secondaryReserve(:, :, 1) = V0;
tempSecondaryReserve(:, :, 1) = V0;
expectedSecondaryReserve(:, :, 1) = V0;
assets(:, :, 1) = V0;

% Setup stress test matrices 5(ii)
secondaryFinCrisisBh100 = zeros(1,12);
secondaryFinCrisisBh100(1) = V0;
assetsFinCrisisBh100 = zeros(1,12);
assetsFinCrisisBh100(1) = V0;
secondaryFinCrisisBh5050 = zeros(1,12);
secondaryFinCrisisBh5050(1) = V0;
assetsFinCrisisBh5050 = zeros(1,12);
assetsFinCrisisBh5050(1) = V0;

%% Reserves and assets
disp('Calculating second order reserves and assets...')
insolvent = zeros(nStrategies, scenarios);

portfolioReturns = zeros(nStrategies, periods - 1, scenarios);
portfolioReturns(1, :, :) = portfolioReturnA(:, :);
portfolioReturns(2, :, :) = portfolioReturnB(:, :);
portfolioReturns(3, :, :) = portfolioReturnC(:, :);
portfolioReturns(4, :, :) = portfolioReturnD(:, :);
portfolioReturns(5, :, :) = portfolioReturnE(:, :);
portfolioReturns(6, :, :) = portfolioReturnF(:, :);

for i = 1:scenarios
    for j = 2:periods
        
        % 5(ii)
        if i == 1 && j <= 12
            % Secondary reserve for buy and hold 100% 15 year bond
            secondaryFinCrisisBh100(j) = calculateSecondaryReserve(secondaryFinCrisisBh100(j-1), ...
                portfolioReturnFinCrisisA(j), P, dG, data(j, 2), data(j, 4));
            
            % Assets for buy and hold 100% 15 year bond
            assetsFinCrisisBh100(j) = calculateAssets(population(j, 9), population(j, 8), population(j, 10), ...
                assetsFinCrisisBh100(i, j-1), portfolioReturnFinCrisisA(j), firstReserve(j-1), secondaryFinCrisisBh100(j-1), ...
                    P, 0);
            
            % Secondary reserve for buy and hold 50% 5 year bond 50% stocks
            secondaryFinCrisisBh5050(j) = calculateSecondaryReserve(secondaryFinCrisisBh5050(j-1), ...
                portfolioReturnFinCrisisB(j), P, dG, data(j, 2), data(j, 4));
            
            % Assets for buy and hold 50% 5 year bond 50% stocks
            assetsFinCrisisBh5050(j) = calculateAssets(population(j, 9), population(j, 8), population(j, 10), ...
                assetsFinCrisisBh5050(i, j-1), portfolioReturnFinCrisisB(j), firstReserve(j-1), secondaryFinCrisisBh5050(j-1), ...
                    P, 0);
        end
       
        % Calculate V, E[V] and assets for all strategies
        for k = 1:nStrategies
            
            % Secondary reserve and payout
            [secondaryReserve(k, i, j), payOut(k, i, j)] = calculateSecondaryReserve(secondaryReserve(k, i, j-1), ...
                portfolioReturns(k, j-1, i), P, dG, data(j, 2), data(j, 4));
            
            % Expected secondary reserve
            expectedSecondaryReserve(k, i, j) = secondaryReserve(k, i, j) * (data(j, 9) + data(j, 10));
            
            % Assets
            assets(k, i, j) = calculateAssets(population(j, 9), population(j, 8), population(j, 10), ...
                assets(k, i, j-1), portfolioReturns(k, j-1, i), firstReserve(j-1), secondaryReserve(k, i, j-1), ...
                    P, payOut(k, i, j));
            if assets(k, i, j) < BE(i, j)
                insolvent(k, i) = 1;
            end
        end
    end
end

%% Calculate insolvency ratios
for i = 1:nStrategies
    nrOfTimesInsolvent(i) = sum(insolvent(i, :) == 1);
    insolvencyRatio(i) = nrOfTimesInsolvent(i) / scenarios; 
end

%% BoF calculation
disp('Calculating BoF and SCR...')

for i = 1:nStrategies
    for j = 1:periods
        bOF(i, :, j) = assets(i, :, j) - BE(:, j)';
    end
end

% SCR strategy a)
deltaBofA = bOF(1, :, 13) - bOF(1, :, 1);
SCRA = quantile(deltaBofA, 0.05);

% SCR strategy b)
deltaBofB = bOF(2, :, 13) - bOF(2, :, 1);
SCRB = quantile(deltaBofB, 0.05);

% SCR strategy c)
deltaBofC = bOF(3, :, 13) - bOF(3, :, 1);
SCRC = quantile(deltaBofC, 0.05);

% SCR strategy d)
deltaBofD = bOF(4, :, 13) - bOF(4, :, 1);
SCRD = quantile(deltaBofD, 0.05);

% SCR strategy e)
deltaBofE = bOF(5, :, 13) - bOF(5, :, 1);
SCRE = quantile(deltaBofE, 0.05);

% SCR strategy f)
deltaBofF = bOF(6, :, 13) - bOF(6, :, 1);
SCRF = quantile(deltaBofF, 0.05);


%% Mean stuff
% bOF5aYr(find(((bOF5aYr==inf)+(bOF5aYr==-inf)+isnan(bOF5aYr))==1)) = 0;
% bOF5bYr(find(((bOF5bYr==inf)+(bOF5bYr==-inf)+isnan(bOF5bYr))==1)) = 0;
% bOF5cYr(find(((bOF5cYr==inf)+(bOF5cYr==-inf)+isnan(bOF5cYr))==1)) = 0;
% bOF5dYr(find(((bOF5dYr==inf)+(bOF5dYr==-inf)+isnan(bOF5dYr))==1)) = 0;
% bOF5eYr(find(((bOF5eYr==inf)+(bOF5eYr==-inf)+isnan(bOF5eYr))==1)) = 0;
% bOF5fYr(find(((bOF5fYr==inf)+(bOF5fYr==-inf)+isnan(bOF5fYr))==1)) = 0;
% 
% meanA = mean(mean(insolventA .* bOF5aYr(:,end)));
% meanB = mean(mean(insolventB .* bOF5bYr(:,end)));
% meanC = mean(mean(insolventC .* bOF5cYr(:,end)));
% meanD = mean(mean(insolventD .* bOF5dYr(:,end)));
% meanE = mean(mean(insolventE .* bOF5eYr(:,end)));
% meanF = mean(mean(insolventF .* bOF5fYr(:,end)));
% 
% %% Plot
% figure
% scatter(insolvencyRatioA, meanA)
% hold on
% scatter(insolvencyRatioB, meanB)
% scatter(insolvencyRatioC, meanC)
% scatter(insolvencyRatioD, meanD)
% scatter(insolvencyRatioE, meanE)
% scatter(insolvencyRatioF, meanF)
% title('Mean-insolvency probability diagram');
% ylabel('Mean of own funds');
% xlabel('Insolvency probability');
% legend('BH(100/0','BH(50/50)', 'CM(50/50)', 'CH(50/50) with duration matching', 'CPPI(1)', 'The Dart throwing strategy');

%% Plot more
plotStuff(expectedSecondaryReserve(1,1,:), assets(1,1,:), BE, bOF(1,1, :), firstReserve, 'Strategy BH(100/0) - scenario 1');
plotStuff(expectedSecondaryReserve(1,2,:), assets(1,2,:), BE, bOF(1,2, :), firstReserve, 'Strategy BH(100/0) - scenario 2');
plotStuff(expectedSecondaryReserve(2,1,:), assets(2,1,:), BE, bOF(2,1, :), firstReserve, 'Strategy BH(50/50) - scenario 1');
plotStuff(expectedSecondaryReserve(2,2,:), assets(2,2,:), BE, bOF(2,2, :), firstReserve, 'Strategy BH(50/50) - scenario 2');
plotStuff(expectedSecondaryReserve(3,1,:), assets(3,1,:), BE, bOF(3,1, :), firstReserve, 'Strategy CM(50/50) - scenario 1');
plotStuff(expectedSecondaryReserve(3,2,:), assets(3,2,:), BE, bOF(3,2, :), firstReserve, 'Strategy CM(50/50) - scenario 2');
plotStuff(expectedSecondaryReserve(4,1,:), assets(4,1,:), BE, bOF(4,1, :), firstReserve, 'Strategy CM(50/50) with duration matching - scenario 1');
plotStuff(expectedSecondaryReserve(4,2,:), assets(4,2,:), BE, bOF(4,2, :), firstReserve, 'Strategy CM(50/50) with duration matching - scenario 2');
plotStuff(expectedSecondaryReserve(5,1,:), assets(5,1,:), BE, bOF(5,1, :), firstReserve, 'Strategy CPPI(1) with duration matching - scenario 1');
plotStuff(expectedSecondaryReserve(5,2,:), assets(5,2,:), BE, bOF(5,2, :), firstReserve, 'Strategy CPPI(1) with duration matching - scenario 2');
plotStuff(expectedSecondaryReserve(6,1,:), assets(6,1,:), BE, bOF(6,1, :), firstReserve, 'Dart Throwing Strategy - scenario 1');
plotStuff(expectedSecondaryReserve(6,2,:), assets(6,2,:), BE, bOF(6,2, :), firstReserve, 'Dart Throwing Strategy - scenario 2');
disp('Finished!')

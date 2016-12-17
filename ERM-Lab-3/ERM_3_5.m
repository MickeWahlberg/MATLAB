%%
% tRange = 60*12:120*12;
% XPVU = zeros(length(tRange),1);
% XPVP = zeros(length(tRange),1);
% i = 1;
% for t = 60*12:120*12;
%     XPVU(i) = PVU(t,dG);
%     XPVP(i) = PVP(t,dG);
%     i = i + 1;
% end

%% Comments
% riskpremie + premie -> avkastning -> utbetalning
% ratios: 5%, 25%, 24% 10% 1% ??
% scr: 38k, 250k, 240k, ....
%% Clear command window
clc
load('data');

%If data not already loaded:
if not(exist('population', 'var'))
    %Load data
    disp('Loading population...')
    f = load('population');
    population = f.population;
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
firstReserve = zeros(721,1);
i = 1;
for t = 60*12:120*12
    firstReserve(i) = G(t, dG);
    i = i + 1;
end
expectedFirstReserve = firstReserve .* (population(:, 9) + population(:, 10));

%% Transform ScenarioFile
adjScenarioData = scenarioData;
adjScenarioData(:,13:17) = exp(-scenarioData(:,13:17).* T);

%% Create returns
a = adjScenarioData(:,3:end);
r1 = 1:length(a);
r2 = r1;
% Remove first month from each scenario
r1(1:121:end) = [];
% Remove last month from each scenario
r2(121:121:end) = [];
% Calculate returns (= R(t)/R(t-1))
normalReturns =  a(r1,:) ./ a(r2,:);

%% Nelson-Siegel and Best Estimate
disp('Creating yield curve...')
% Get a 60-year yield curve with nelson-siegel
rate = zeros(121,2000, 720);
for t =  1:121
    for i = 1:2000
        beta = calibrateNSParameters(T, scenarioData(t+((i-1)*121) , 13:17)');
        rate(t, i, 1:end-(t-1)) = nelson(beta, 1/12 : 1/12 : ((60*12-t+1)/12));
    end
end

%%
% First order cash flow
firstOrderCashFlow = [0; b .* population(2:end, 10) ./ 12 + firstReserve(1:end-1) .* ...
    population(2:end, 8) - P .* population(2:end, 9) ./ 12]';

% Calculate BE
BE = zeros(2000, 121);
for i = 1:121
   t = 1:(721-i);
   rates = squeeze(rate(i,:,t));
   BE(:,i) = sum(firstOrderCashFlow((i+1):end).*exp(-rates.*(t./12)),2);
end
%%
% Calculate ZCB-investing with "rolling over"
r = scenarioData(:, 13:end);
t1 = 1:(scenarios * periods);
t0 = t1;

% Remove first month from each scenario
t1(1:121:end) = [];
% Remove last month from each scenario
t0(121:121:end) = [];

zcbReturns = exp(-(r(t1,:) - r(t0,:)) .* T) .* exp(r(t1,:) * 1/12);
normalReturns(:, 11:end) = zcbReturns;

%% Calculate Returns
nStrategies = 6;
portfolioReturns = zeros(nStrategies, periods - 1, scenarios);
for i = 1:nStrategies
    switch i
        case 1
            %% Assignment A
            portfolioReturns(i, :, :) = calculatePortfolioReturnA(normalReturns);
        case 2
            %% Assignment B
            portfolioReturns(i, :, :) = calculatePortfolioReturnB(normalReturns);
        case 3
            %% Assignment C
            disp('Calculating portfolio returns with CM(50/50)-strategy..')
            portfolioReturns(i, :, :) = calculatePortfolioReturnC(normalReturns);
        case 4
            %% Assignment D:
            %CM(50/50) with bond investments matching liability duration in t=0.
            disp('Calculating portfolio returns with CM(50/50) duration matching strategy..')
            portfolioReturns(i, :, :) = calculatePortfolioReturnD(adjScenarioData, rate, BE, firstOrderCashFlow, normalReturns);
        case 5
            %% Assignment E:
            portfolioReturns(i, :, :) = calculatePortfolioReturnE(adjScenarioData, rate, BE, firstOrderCashFlow, normalReturns);
        case 6
            %% Assignment F:
            %The Dart Strategy
            portfolioReturns(i, :, :) = calculatePortfolioReturnF(normalReturns);
    end
end
%% Create returns for stress test
[portfolioReturnFinCrisisA, unAdjFinancialCrisisData] = calculateFinancialCrisisReturn(marketData, 'BH100');
portfolioReturnFinCrisisB = calculateFinancialCrisisReturn(marketData, 'BH5050');

%% ALM data

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
finCrisisBE = zeros(1,12);
rejt = zeros(12,720);
%% Reserves and assets
disp('Calculating second order reserves and assets...')
insolvent = zeros(nStrategies, scenarios);

for i = 1:scenarios
    for j = 2:periods
        
        % 5(ii)
        if i == 1 && j <=13
            t = 1:(722-j);
            betaFinCrisis = calibrateNSParameters(T, unAdjFinancialCrisisData(j, 11:15)');
            rateFinCrisis = nelson(betaFinCrisis, 1/12 : 1/12 : 60);
            finCrisisBE(j-1) = sum(firstOrderCashFlow((j):end).*exp(-rateFinCrisis(1:end-(j-2)).*(t./12)));
        end
        
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

        end
    end
end


%% Insolvency
insolvencyRatio = zeros(1,nStrategies);
for i = 1:nStrategies
   insolvencyRatio(i) = mean(sum(BE>squeeze(assets(i,:,:)),2)>0);
end

%% BoF calculation
disp('Calculating BoF and SCR...')

for i = 1:nStrategies
        bOF(i, :, :) = squeeze(assets(i, :, :)) - BE(:, :);
end
deltaBof = zeros(nStrategies, scenarios);
SCR = zeros(1, nStrategies);

alpha = 0.05;

for i = 1:nStrategies
   deltaBof(i,:) = bOF(i, :, 13) - bOF(i, :, 1);
   SCR(i) = quantile(deltaBof(i,:), alpha);
end

finCrisisBh100DeltaBOF = (assetsFinCrisisBh100(end) - finCrisisBE(end)) - (assetsFinCrisisBh100(1) - finCrisisBE(1));
finCrisisBh5050DeltaBOF = (assetsFinCrisisBh5050(end) - finCrisisBE(end)) - (assetsFinCrisisBh5050(1) - finCrisisBE(1));
%% Mean stuff

figure
hold on
for i = 1:nStrategies
    notInsolvent = sum(BE<squeeze(assets(i,:,:)),2)==121;
    Y = squeeze(bOF(i,:,end));
    Y = mean(Y(1,notInsolvent));
    
    scatter(insolvencyRatio(i), Y, 'filled')
end
title('Mean-insolvency probability diagram');
ylabel('Mean of own funds');
xlabel('Insolvency probability');
legend('BH(100/0','BH(50/50)', 'CM(50/50)', 'CH(50/50) with duration matching', 'CPPI(1)', 'The Dart throwing strategy');
hold off



%% Plot more
titles = {'BH(100/0)', 'BH(50/50)', 'CM(50/50)', 'CM(50/50) with duration matching', 'CPPI(1) with duration matching', 'Dart Throwing'};
for strategy = 1:nStrategies
   for i = 1:2
      titleString =  sprintf('Strategy %s - scenario %d', titles{strategy}, i);
      plotStuff(expectedSecondaryReserve(strategy,i,:), assets(strategy,i,:), BE(i,:), bOF(strategy,i, :), expectedFirstReserve(1:121), titleString);
   end
end
disp('Finished!')

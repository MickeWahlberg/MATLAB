%% Load data
%If data not already loaded:
if not(exist('scenarioData', 'var'))
    %Load data
    disp('Loading scenario data...')
    %f = load('ScenarioFile_dcc_labb3');
    scenarioData = xlsread('ScenarioFile_Excel');
    scenarioData = scenarioData(1:121*2000,:);
    disp('Scenariod data loaded!')
end

%% Efficient frontiers of three risk measures:
% Variance, LPM2, and ES95%

% Calculate stock returns
stocks = scenarioData(:,3:12);
R = stocks(13:121:end,:) ./ stocks(1:121:end,:);

% Number of scenarios and stocks
[nTimeSteps, nStocks] = size(stocks);
nScenarios = nTimeSteps/121;
% Calculate covariance matrix (sigma) and mean return (Mu = m)
Rstar = R - mean(R);
sigma = (1 / nScenarios) * (Rstar' * Rstar);
m = mean(R)';
stockVol = sqrt(diag(sigma));

%Starting Capital
V0 = 1;

s = length(sigma);
I = ones(s,1);

%% LPM parameters for quadprog
% Length of X (which quadprog solves for) and tau (= starting capital, V0)
L1 = nScenarios + nStocks;
tau = V0;
% A (for A*X <= b)
ALPM = -1 * eye(L1);
ALPM(1:nScenarios, (end - (nStocks - 1)):end) = -1 * R;
ALPM((nScenarios + 1):end, :) = 0;
% b
bLPM = -tau * ones(1, L1);
bLPM((nScenarios + 1):end) = 0;
% H (for X'*H*X)
H = (1 / nScenarios) * eye(L1);
H((nScenarios + 1):end, :) = 0;
% m (Mu) 
mLPM = zeros(L1, 1);
mLPM((nScenarios + 1):end) = m;
% I (array of ones, and in this case zeroes, for summation of weights)
ILPM = ones(L1, 1);
ILPM(1:nScenarios) = 0;

%% ES parameters for quadprog
alpha = 0.95;
% Length of X (which quadprog solves for) and tau (= starting capital, V0)
L2 = nScenarios + nStocks + 1;
% A (for A*X <= b)
AES = -1 * eye(L2);
AES(:,nScenarios + 1) = -1;
AES(1:nScenarios, (nScenarios + 2):end) = -R + 1;
AES(L2:end, :) = 0;
% b
bES = zeros(1, L2);
% H (for X'*H*X)
f = (1 / ((1 - alpha) * nScenarios)) * ones(L2,1);
f(nScenarios + 1) = 1;
f((nScenarios + 2):end) = 0;
% m (Mu) 
mES = zeros(L2, 1);
mES((nScenarios + 2):end) = m;
% I (array of ones, and in this case zeroes, for summation of weights)
IES = ones(L2, 1);
IES(1:(nScenarios + 1)) = 0;





%%
% Values of m0 (Mu Zero) to test
m0Range = 1.09:0.005:1.17;

% Number of simulations
l = length(m0Range);

% Variables for holding weights and volatilities of risk measures
% Variance
wVar = zeros(10,l);
volVar = zeros(l,1);
% LPM
XLPM = zeros(L1,l);
wLPM = zeros(10,l);
volLPM = zeros(l,1);
% Expected Shortfall
XES = zeros(L2,l);
wES = zeros(10,l);
volES = zeros(l,1);

index = 1;
 
% Option for using the dual-simplex method in linprog
options = optimoptions('linprog','Algorithm','dual-simplex');

% For all m0:s
for m0 = m0Range
    
    % Variance
    wVar(:,index) = quadprog(sigma, [], [], [], [I';m'], [V0;m0*V0], zeros(1,nStocks));
    volVar(index) = sqrt(wVar(:,index)' * sigma * wVar(:,index));
  
    % LPM (Beta = 2) 
    XLPM(:,index) = quadprog(H, [], ALPM, bLPM', [ILPM';mLPM'], [V0;m0*V0], zeros(1,L1));
    wLPM(:,index) = XLPM((nScenarios + 1):end,index);
    volLPM(index) = sqrt(wLPM(:,index)' * sigma * wLPM(:,index));
    
    % ES - Expected Shortfall
    XES(:,index) = linprog(f, AES, bES', [IES';mES'], [V0;m0*V0], zeros(1,L2), [], options);
    wES(:,index) = XES((nScenarios + 2):end,index);
    volES(index) = sqrt(wES(:,index)' * sigma * wES(:,index));
    
    index = index + 1;
end

%% Plot
% Mean-volatility plot with individual stocks and efficient frontier
Y = m0Range;
figure
plot(stockVol, m, 'o');
hold on 
plot(volVar,Y);
plot(volLPM, Y);
plot(volES, Y);
title('Mean-volatility Plot With Individual Stocks and Efficient Frontier');
legend('Stocks','Variance', 'LPM', 'ES');
hold off

% Area plots of stock weights for risk measures 
figure

subplot(2,2,1);
area(wVar')
title('Area Plot for Variance as Risk Measure');
subplot(2,2,2);
area(wLPM')
title('Area Plot for LPM (Beta = 2) as Risk Measure');
subplot(2,2,3);
area(wES')
title('Area Plot for Expected Shortfall as Risk Measure');
legend('Boeing', 'Caterpillar', 'Coca Cola', 'Dupont', 'JP Morgan', '3m', 'Microsoft', 'Pfizer', 'Walmart', 'Exxon Mobile');
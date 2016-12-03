%%Load data
load('MarketData9615');

%% Assignment 1
stocks = MarketData(:,1:10);
V0 = 1;

% Returns
y = stocks(2:end,:)./stocks(1:(end-1),:) - 1 ;

m = mean(y(1:100,:));

% Portfolio covariance, variance, and volatility
portfolioCov = cov(y(1:100,:));
portfolioVariance = portfolioCov(1:11:end);
portfolioVolatility = sqrt(portfolioVariance);

%Estimate parameters
s = length(portfolioCov);
I = ones(s,1);

% For same dimensions as in notes
m = m';

A = I' / portfolioCov * I;
B = I' / portfolioCov * m;
C = m' / portfolioCov * m;
D = A * C - B^2;

% Efficient frontier
% Short selling allowed
M = -0.02:0.001:0.03;
efficientFrontierI = sqrt((A/D) .* M.^2 - (2 * B * V0 .* M) / D + C * V0^2 / D);

% Short selling NOT allowed

MII = 0.001:0.001:0.03;
index = 1;
weightsEfficII = zeros(length(MII),10);
efficientFrontierII = zeros(1,length(MII));

for i = MII
    weightsEfficII(index,:) = quadprog(portfolioCov, [], [], [], [ones(1,s);m'], [V0;i], zeros(1,s));
    efficientFrontierII(index) = sqrt(weightsEfficII(index, :) * portfolioCov * weightsEfficII(index, :)');
    index = index + 1;
end

%% Sharpe, short selling allowed
sharpeAlloc= V0 * (m' / portfolioCov) / ((I' / portfolioCov) * m);
sharp = sqrt(sharpeAlloc * portfolioCov * sharpeAlloc');
m_sharp = sharpeAlloc * m;

% %% Sharpe, short selling NOT allowed
minVarAlloc = ((A * mean(m) * V0 - B * V0) / D) * (portfolioCov \ m) + ((C * V0 - B * mean(m) * V0) / D) * (portfolioCov \ I);
minVarPortfolio = sqrt(minVarAlloc' * portfolioCov * minVarAlloc);
m_minVar = minVarAlloc' * m;

%% Plots

plot(portfolioVolatility,m,'o');
hold on 
plot(efficientFrontierI, M);
plot(efficientFrontierII, MII);
plot(sharp, m_sharp, 'x');


% Plot tangent
X = linspace(0.0,0.07);
y1 = m_sharp;
Y = ((y1-0)*X+0*sharp-y1*0)/(sharp-0);
plot(X,Y)
hold off


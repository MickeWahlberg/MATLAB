
%%Load data
if not(exist('sp500Struct', 'var'))
    %Load data
    disp('Loading SP500...')
    sp500Struct = hist_stock_data('01011996','31122015', '^gspc','frequency', 'w');
    sp500 = sp500Struct.AdjClose;
    disp('Finished!')  
end

if not(exist('MarketData', 'var'))
    %Load data
    disp('Loading Stock Data...')
    load('MarketData9615');
    stocks = MarketData(:,1:10);
    disp('Finished!')  
end

%% Benchmark allocation



% Return on stocks + risk-free asset (R0 = 1) and benchmark (S&P 500)
yStocks = stocks(2:end,:)./stocks(1:(end-1),:);
yAssets = yStocks;
ySP500 = sp500(2:end,:)./sp500(1:(end-1),:);

% Timesteps

% Starting fund in assets and benchmark
V0 = 1;
L0 = V0;
R0 = 1;


% Fetch 100 last returns of stocks and benchmark (S&P 500)
yAssets = yAssets((end-99):end,:);
ySP500 = ySP500((end-99):end,:);


m = mean(yAssets)';
I = ones(size(m));

% Returns on assets and benchmark
R = [yAssets, ySP500];

%R de-meaned
Rdm = R - mean(R);

% MA-estimation of covariance matrix of R
T = 100;
epsilonR = (1/T) * (Rdm' * Rdm);

%Covariance of assets and benchmark (L)
sigma = epsilonR(1:(end-1),1:(end-1));
sigmaLR = L0 * epsilonR(1:(end-1),end);
wh = sigma \ sigmaLR;


% Range of expected return
m0Range = 1.00001:0.0001:(max(m)+0.01);


s = length(m0Range);
column = 1;

% For holding asset weights
wstar = zeros(10, s);
w = zeros(10, s);
w0 = zeros(1, s);
volatilities = zeros(s,1);

     
for m0 = m0Range
    % Calclulate weights of stocks and risk-free asset
    lambda1 = (wh' * (m - R0 * I) + V0 * (R0  - m0)) / ((m - R0 * I)' * (sigma \ (m - R0 * I)));
    w(:,column) = -lambda1 * (sigma \ (m - R0 *I))   + wh;
    w0(column) = V0 - w(:,column)' * I;
    
    % Calculate volatility
    volatilities(column) = sqrt(w(:,column)' * sigma * w(:,column));
    
    column = column + 1;
end
   
RV = [w' .* R(end,1:(end-1)), w0' * R0];
V1 = sum(RV,2);
L = L0 * R(end,end);
meanTrackingError = V1 - L;

%% Plots
portfolioVariance = diag(sigma);
portfolioVolatility = sqrt(portfolioVariance);

figure
plot(portfolioVolatility,m,'o');
title('Mean-volatility chart of stocks and efficient frontier');
%legend('Boeing', 'Caterpillar', 'Coca Cola', 'Dupont', 'JP Morgan', '3m', 'Microsoft', 'Pfizer', 'Walmart', 'Exxon Mobile');
hold on
plot(volatilities,m0Range);
hold off

figure 
plot(meanTrackingError,m0Range);


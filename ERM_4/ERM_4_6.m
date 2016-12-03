%% Load data
%If data not already loaded:
if not(exist('scenarioData', 'var'))
    %Load data
    disp('Loading scenario data...')
    f = load('scenarioData');
    scenarioData = f.scenarioData;
    disp('Scenariod data loaded!')
end

%% Efficient frontiers of three risk measures:
% Variance, LPM2, and ES95%

%Calculate stock returns
stocks = scenarioData(:,3:12);
[rows, columns] = size(stocks);
y = zeros(rows - 1 * 2000, columns);
nScenarios = 2000;
yMean = zeros(nScenarios * 10,10);
sigma = zeros(10);
for i = 1:nScenarios
    startIndex = 1 + 121 * (i - 1);
    endIndex = 121 * i;
    s = stocks(startIndex:12:endIndex, :);
    y = s(2:end,:) ./ s(1:(end-1),:);
    sigma = sigma +  cov(y);
    yMean(i,:) = mean(y);
end

m = mean(yMean);
sigma = sigma ./ nScenarios;
V0 = 1;




m = mean(y)';
s = length(sigma);
I = ones(s,1);




l = length(1.09:0.005:1.17);
wVar = zeros(10,l);
index = 1;
volVar = zeros(l,1);
for m0 = 1.09:0.005:1.17
    
    if min(m) > m0
        wVar(:, index) = (m == min(m)).*V0;
    elseif max(m) < m0
        wVar(:, index) = (m == max(m)).*V0;
    else
        wVar(:,index) = quadprog(sigma, [], [], [], [I';m'], [V0;m0], zeros(1,s));
    end
    
    volVar(index) = sqrt(wVar(:,index)' * sigma * wVar(:,index));
    
    
    
    
    
    index = index + 1;
end

Y = 1.09:0.005:1.17;
plot(volVar,Y);

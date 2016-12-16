function portfolioReturn = calculatePortfolioReturnC(normalReturns)
wTest = [ones(1,10)/20 0.5];
portfolioReturn = sum(normalReturns(:,1:11) .*  wTest,2) - 1;
portfolioReturn = reshape(portfolioReturn,120, length(portfolioReturn)/120);
end
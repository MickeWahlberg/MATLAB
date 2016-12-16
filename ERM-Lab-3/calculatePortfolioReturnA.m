function portfolioReturn = calculatePortfolioReturnA(normalReturns)
portfolioReturn = normalReturns(:,13) - 1;
portfolioReturn = reshape(portfolioReturn,120, length(portfolioReturn)/120);
end
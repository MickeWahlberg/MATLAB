function portfolioReturn = calculatePortfolioReturnF(normalReturns)
%The total capital is invested in a random asset in each time step.
portfolioReturnF = zeros(120, 2000);
asset = randi([1 15], 120,1);
for t = 1:120
    portfolioReturnF(t,:) = normalReturns(t:120:end, asset(t));
end

portfolioReturn = portfolioReturnF - 1;
end
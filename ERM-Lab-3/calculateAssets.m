function assets = calculateAssets(premiumPayers, transfers, retired, assets, portfolioReturn, firstReserve, secondaryReserve, P, payOut)
        %  Calculate value of assets
        assetReturns = portfolioReturn * assets;
        vCashFlows = (P / 12) * premiumPayers - transfers * ...
            max(firstReserve, secondaryReserve) - payOut * retired;
        assets = assets + vCashFlows + assetReturns;
end

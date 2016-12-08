function assets = calculateAssets(premiumPayers, transfers, retierd, assets, portfolioReturn, firstReserve, secondaryReserve, P, payOut)
        %  Calculate value of assets
        assetReturns = (1 + portfolioReturn) * assets - assets;
        vCashFlows = (P / 12) * premiumPayers - transfers * ...
            max(firstReserve, secondaryReserve) - ...
                payOut * retierd;
        assets = assets + vCashFlows + assetReturns;
end

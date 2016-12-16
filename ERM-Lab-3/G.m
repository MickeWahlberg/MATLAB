%First order reserve for year x
function firstOrderReserve = G(x, dG)
    b = 170619.01695;
    P = 120000;
    firstOrderReserve = (b / 12) .* PVU(x, dG) - (P / 12) .* PVP(x,dG);
end
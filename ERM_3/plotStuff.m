function plota = plotStuff(expectedSecondaryReserve, assets, BE, bOF, firstReserve, titlen)
tempExp(:,:) = expectedSecondaryReserve(1,1,:);
tempAss(:,:) = assets(1, 1, :);
tempBOF(:,:) = bOF(1, 1, :);
figure
plot(tempExp);
hold on
plot(tempAss);
plot(BE(1, 1:121));
plot(tempBOF);
plot(firstReserve(1:121));
title(sprintf(titlen));
legend('Secondary reserve', 'Assets', 'Best estimate', 'BOF', 'First reserve');
plota = 1;
end
function result = PVU(x, dG)
    %Retirement and time to retirement
    z = 65*12;
    u = z - x;
    S = 12*(120 - 65);

    t = (u*(u>0) + 1):(u + S);
    result = sum((L((x + t)./12) ./ L(x/12)) .* exp(-dG .* t/12));
end
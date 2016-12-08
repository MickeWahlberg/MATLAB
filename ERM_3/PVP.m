%Returns present value of premiums (1 kr)
function result = PVP(x, dG)
    %Retirement and time to retirement
    z = 65*12;
    u = z - x;
    
        t = 1: u;
        result = (u>0).*sum((L((x + t)./12) ./ L(x/12)) .* exp(-dG .* t/12));
end
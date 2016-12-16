function result = PVUII(x, dG)
    %Retirement and time to retirement
    z = 65;
    u = z - x;
    S = 120 - 65;
   
    
    if x > 65
        t = 0:(1/12): u + S;
        result = sum((L(x + t) ./ L(x)) .* exp(-dG * t));
    else 
        result = 0;
    end
      
end
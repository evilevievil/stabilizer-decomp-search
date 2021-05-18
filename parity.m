function par = parity(x)
    par = uint8(0);
    for i = 1:8
        par = par + bitget(x,i);
    end
    par = bitget(par,1);
end
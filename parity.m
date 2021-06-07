function par = parity(x)
    par = uint8(0);
    for i = 1:8
        par = bitxor(par,bitget(x,i));
    end
end

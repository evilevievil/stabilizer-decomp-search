function par = parity(x)
    par = const.init_uint;
    for i = 1:16
        par = bitxor(par,bitget(x,i));
    end
end

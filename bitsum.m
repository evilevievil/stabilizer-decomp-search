function bit_sum = bitsum(x)
    bit_sum = const.init_uint;
    for i = 1:16
        bit_sum = bit_sum + bitget(x,i);
    end
end

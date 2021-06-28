function bit_sum = bitsum(x)
    bit_sum = const.init_uint;
    for i = 1:8
        bit_sum = bit_sum + bitget(x,i);
    end
end

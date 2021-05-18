function bit_sum = bitsum(x)
    bit_sum = uint8(0);
    for i = 1:8
        bit_sum = bit_sum + bitget(x,i);
    end
end